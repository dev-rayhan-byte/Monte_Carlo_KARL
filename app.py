import streamlit as st
import pandas as pd
import numpy as np
import threading
import time
import os
import math
import traceback
import imageio.v2 as imageio
from openpyxl import Workbook
from ase.cluster import FaceCenteredCubic, BodyCenteredCubic, HexagonalClosedPacked
from ase.io import write, read
from ase.neighborlist import build_neighbor_list
import matplotlib.pyplot as plt
from ase.visualize.plot import plot_atoms
import base64

# --- Page and Session State Configuration ---
st.set_page_config(page_title="Monte Carlo Nanoparticle Simulator", layout="wide")

if 'simulation_results' not in st.session_state:
    st.session_state.simulation_results = None
if 'simulation_running' not in st.session_state:
    st.session_state.simulation_running = False
if 'error_message' not in st.session_state:
    st.session_state.error_message = None

# --- Core Simulation & Visualization Functions (Adapted from your script) ---
BOLTZMANN_K = 8.617333262e-5
BULK_COORD = 12 # Coordination number for bulk atoms in FCC

lattice_map = {
    'fcc': FaceCenteredCubic,
    'bcc': BodyCenteredCubic,
    'hcp': HexagonalClosedPacked
}

# --- Helper Functions ---
def symbol_type(sym, A):
    return 'A' if sym == A else 'B'

def calculate_energy(p, A, coeffs):
    """Calculates the total energy of the nanoparticle."""
    nl = build_neighbor_list(p, bothways=True, self_interaction=False)
    count = {k: 0 for k in coeffs}
    for i in range(len(p)):
        t_i = symbol_type(p[i].symbol, A)
        neighbors = nl.get_neighbors(i)[0]
        if len(neighbors) < BULK_COORD:
            count[f'x{t_i}-S'] += 1
        for j in neighbors:
            if i < j:
                t_j = symbol_type(p[j].symbol, A)
                # Sort to handle A-B vs B-A consistently
                pair = '-'.join(sorted([t_i, t_j]))
                key = f'x{pair}'
                count[key] += 1
    return sum(coeffs[k] * count.get(k, 0) for k in coeffs)


def count_surface(p, A):
    """Counts surface atoms and the ratio of element A on the surface."""
    nl = build_neighbor_list(p, bothways=True, self_interaction=False)
    total_A, surf, surf_A = 0, 0, 0
    for i in range(len(p)):
        neighbors = nl.get_neighbors(i)[0]
        t = symbol_type(p[i].symbol, A)
        if t == 'A':
            total_A += 1
        if len(neighbors) < BULK_COORD:
            surf += 1
            if t == 'A':
                surf_A += 1
    ratio = surf_A / surf if surf else 0
    return total_A, surf, surf_A, ratio

def run_monte_carlo_simulation(params, progress_callback):
    """
    Main function to run the MC simulation.
    Saves structures, log, and trajectory files.
    """
    A = params['element_A']
    B = params['element_B']
    T = params['temperature']
    N_STEPS = params['n_steps']
    SAVE_INTERVAL = params['save_interval']
    SNAPSHOT_INTERVAL = params['snapshot_interval']

    # Generate initial particle
    ClusterBuilder = lattice_map.get(params['lattice_type'])
    particle = ClusterBuilder(A, surfaces=params['surfaces'], layers=params['layers'])
    n_atoms = len(particle)
    n_A = int(n_atoms * params['composition_A'])
    indices_A = np.random.choice(range(n_atoms), size=n_A, replace=False)
    for i in range(n_atoms):
        particle[i].symbol = A if i in indices_A else B

    # Create directories for output
    os.makedirs("trajectory", exist_ok=True)
    os.makedirs("output", exist_ok=True)
    
    initial_xyz_path = f"output/initial_{A}{B}_{n_atoms}.xyz"
    write(initial_xyz_path, particle)

    # Simulation loop
    energy = calculate_energy(particle, A, params['coefficients'])
    log = []
    start_time = time.time()

    for step in range(1, N_STEPS + 1):
        i = np.random.randint(0, n_atoms)
        neighbors = build_neighbor_list(particle).get_neighbors(i)[0]
        if not neighbors.any(): continue
        j = np.random.choice(neighbors)
        if particle[i].symbol == particle[j].symbol: continue

        trial = particle.copy()
        trial[i].symbol, trial[j].symbol = trial[j].symbol, trial[i].symbol
        dE = calculate_energy(trial, A, params['coefficients']) - energy

        if dE < 0 or np.random.random() < math.exp(-dE / (BOLTZMANN_K * T)):
            particle = trial
            energy += dE

        if step % SAVE_INTERVAL == 0:
            total_A, surf, surf_A, ratio = count_surface(particle, A)
            log_entry = {
                'Step': step, 'Energy (eV)': energy, f'Total {A}': total_A,
                'Surface Atoms': surf, f'{A} on Surface': surf_A, f'Surface {A} Ratio': ratio
            }
            log.append(log_entry)
            progress_callback(step / N_STEPS, f"Step: {step}/{N_STEPS} | Energy: {energy:.3f} eV")
            
        if step % SNAPSHOT_INTERVAL == 0:
            write(f"trajectory/step_{step:05d}.xyz", particle)

    # Save final results
    final_xyz_path = f"output/final_{A}{B}_{n_atoms}.xyz"
    write(final_xyz_path, particle)

    xlsx_path = f"output/MMC_{A}{B}_log.xlsx"
    wb = Workbook()
    ws = wb.active
    ws.title = "Simulation Log"
    meta = [(k, v) for k, v in params.items() if k != 'coefficients']
    for i, (k, v) in enumerate(meta, 1):
        ws.cell(row=i, column=1, value=str(k))
        ws.cell(row=i, column=2, value=str(v))
    
    if log:
        header_row = len(meta) + 2
        headers = list(log[0].keys())
        for j, h in enumerate(headers, 1):
            ws.cell(row=header_row, column=j, value=h)
        for i, entry in enumerate(log, header_row + 1):
            for j, h in enumerate(headers, 1):
                ws.cell(row=i, column=j, value=entry[h])
    wb.save(xlsx_path)

    return {
        "initial_xyz": initial_xyz_path,
        "final_xyz": final_xyz_path,
        "log": pd.DataFrame(log),
        "xlsx_file": xlsx_path,
        "duration": time.time() - start_time,
        "initial_surface_data": count_surface(read(initial_xyz_path), A),
        "final_surface_data": count_surface(particle, A),
        "params": params,
        "particle": particle
    }

def create_evolution_video(log_df, trajectory_folder, output_file, params, progress_callback):
    """Generates an MP4 movie from the simulation trajectory."""
    xyz_files = sorted([f for f in os.listdir(trajectory_folder) if f.endswith(".xyz")])
    if not xyz_files:
        st.warning("No trajectory files found to create a video.")
        return None
        
    images = []
    symbol_to_color = {'Pt': 'silver', 'Ru': 'royalblue'} # Simplified

    for i, xyz_file in enumerate(xyz_files):
        progress_callback(i / len(xyz_files), f"Rendering frame {i+1}/{len(xyz_files)}...")
        atoms = read(os.path.join(trajectory_folder, xyz_file))
        label = xyz_file.replace(".xyz", "")
        current_step = int(''.join(filter(str.isdigit, label)))

        colors = [symbol_to_color.get(atom.symbol, 'gray') for atom in atoms]
        element_counts = {s: atoms.get_chemical_symbols().count(s) for s in np.unique(atoms.get_chemical_symbols())}

        fig, ax = plt.subplots(figsize=(8, 8), dpi=150)
        fig.patch.set_facecolor('#f0f0f0')
        
        plot_atoms(atoms, ax, radii=0.7, colors=colors, rotation=('30x,30y,0z'))
        ax.axis('off')
        
        # Add legend
        x0, y0, dy = 0.05, 0.95, 0.04
        for idx, (symbol, count) in enumerate(sorted(element_counts.items())):
            fig.text(x0, y0 - idx * dy, f"● {symbol} ({count})", fontsize=12, ha='left', va='top',
                     color=symbol_to_color.get(symbol, 'gray'), fontweight='bold')

        # Add inset plot
        inset_ax = fig.add_axes([0.62, 0.72, 0.35, 0.22])
        log_subset = log_df[log_df['Step'] <= current_step]
        if not log_subset.empty:
            inset_ax.plot(log_subset["Step"], log_subset["Energy (eV)"], color='blue', label='Energy')
            ax2 = inset_ax.twinx()
            ax2.plot(log_subset["Step"], log_subset[f"Surface {params['element_A']} Ratio"], color='orange', label='Surface Ratio')
            inset_ax.legend(loc='upper left', fontsize=8)
            ax2.legend(loc='lower left', fontsize=8)
        inset_ax.set_xlabel('MC Step', fontsize=9)
        inset_ax.tick_params(axis='both', which='major', labelsize=8)
        inset_ax.grid(True, linestyle='--', alpha=0.6)

        tmp_png = f"/tmp/{xyz_file}.png"
        fig.savefig(tmp_png)
        plt.close(fig)
        images.append(imageio.imread(tmp_png))

    progress_callback(1.0, "Compiling video...")
    imageio.mimsave(output_file, images, fps=2)
    return output_file

# --- UI Layout ---
st.title("🔬 Monte Carlo Nanoparticle Simulator")

# Sidebar for simulation parameters
st.sidebar.header("1. Simulation Parameters")
with st.sidebar.expander("Elements & Composition", expanded=True):
    element_A = st.text_input("Element A", "Pt")
    element_B = st.text_input("Element B", "Ru")
    composition_A = st.slider("Composition of A", 0.0, 1.0, 0.5)

with st.sidebar.expander("Thermodynamics & Steps"):
    temperature = st.number_input("Temperature (K)", 1.0, 5000.0, 250.0)
    n_steps = st.number_input("MC Steps", 100, 50000, 10000, 100)
    save_interval = st.number_input("Log Interval", 1, 1000, 100)
    snapshot_interval = st.number_input("Snapshot Interval (for video)", 1, 1000, 500)

with st.sidebar.expander("Geometry & Lattice"):
    lattice_type = st.selectbox("Lattice Type", ["fcc", "bcc", "hcp"])
    layers = tuple(st.multiselect("Layers (x,y,z)", list(range(1,11)), [7,7,7]))
    surfaces_str = st.text_input("Surfaces", "[(1,1,1),(1,1,1),(1,1,0)]")

with st.sidebar.expander("Energy Coefficients (eV)"):
    coeffs = {
        'xA-A': st.number_input("A-A Bond", value=-0.022078, format="%.6f"),
        'xB-B': st.number_input("B-B Bond", value=-0.150000, format="%.6f"),
        'xA-B': st.number_input("A-B Bond", value=-0.109575, format="%.6f"),
        'xA-S': st.number_input("A-Surface", value=-0.250717, format="%.6f"),
        'xB-S': st.number_input("B-Surface", value=-0.300000, format="%.6f")
    }

st.sidebar.header("2. Simulation Controls")
run_button = st.sidebar.button("▶️ Run Simulation", type="primary", disabled=st.session_state.simulation_running)
generate_video = st.sidebar.checkbox("Generate Evolution Video", value=True)

# Main area for status and results
status_placeholder = st.empty()
progress_bar = st.sidebar.progress(0, text="Awaiting simulation...")

def simulation_target(params):
    """Target function for the simulation thread."""
    try:
        def progress_cb(pct, text):
            # This callback will be called from the thread
            # We store the progress and let the main thread update the UI
            st.session_state.progress_pct = pct
            st.session_state.progress_text = text
            
        results = run_monte_carlo_simulation(params, progress_cb)
        st.session_state.simulation_results = results
        if generate_video:
            video_progress_bar = st.sidebar.progress(0, text="Starting video generation...")
            def video_progress_cb(pct, text):
                st.session_state.progress_pct = pct
                st.session_state.progress_text = text
                
            video_path = create_evolution_video(
                results['log'], "trajectory", "output/evolution.mp4", params, video_progress_cb
            )
            st.session_state.simulation_results["video_file"] = video_path
    except Exception as e:
        st.session_state.error_message = traceback.format_exc()
    finally:
        st.session_state.simulation_running = False

# --- Button Logic ---
if run_button:
    try:
        surfaces = eval(surfaces_str)
        params = {
            'element_A': element_A, 'element_B': element_B, 'composition_A': composition_A,
            'temperature': temperature, 'n_steps': n_steps, 'save_interval': save_interval,
            'snapshot_interval': snapshot_interval, 'layers': layers, 'surfaces': surfaces,
            'coefficients': coeffs, 'lattice_type': lattice_type
        }
        st.session_state.simulation_running = True
        st.session_state.simulation_results = None
        st.session_state.error_message = None
        st.session_state.progress_pct = 0
        st.session_state.progress_text = "Initializing..."
        
        thread = threading.Thread(target=simulation_target, args=(params,))
        thread.start()
        
    except Exception as e:
        st.error(f"Invalid input parameters: {e}")
        st.session_state.simulation_running = False

# --- Display Loop (updates UI based on session_state) ---
if st.session_state.simulation_running:
    progress_bar.progress(st.session_state.get('progress_pct', 0), text=st.session_state.get('progress_text', ''))
    status_placeholder.info(f"⚙️ {st.session_state.get('progress_text', 'Simulation in progress...')}")
    time.sleep(0.1) # Short sleep to prevent crazy looping
    st.rerun() # Rerun to check thread status and update progress

if st.session_state.error_message:
    st.error("❌ Simulation Failed!")
    with st.expander("Click to see error details"):
        st.code(st.session_state.error_message)

if st.session_state.simulation_results:
    res = st.session_state.simulation_results
    params = res['params']
    A = params['element_A']
    df_log = res['log']

    st.success(f"✅ Simulation Complete! ({res['duration']:.2f} seconds)")

    # --- 1. Summary & Validation ---
    st.subheader("📊 Summary & Validation")
    col1, col2, col3 = st.columns(3)
    init_total, init_surf, init_surf_A, init_ratio = res["initial_surface_data"]
    final_total, final_surf, final_surf_A, final_ratio = res["final_surface_data"]
    col1.metric(f"Initial Surface {A}%", f"{init_ratio:.3%}", help=f"{init_surf_A} of {init_surf} surface atoms")
    col2.metric(f"Final Surface {A}%", f"{final_ratio:.3%}", f"{final_ratio - init_ratio:.3%}", delta_color="inverse")
    col3.metric("Total Atoms", len(res['particle']))
    
    with st.expander("🔍 Validation Checks & Detailed Comparison"):
        st.markdown("---")
        st.markdown("#### **Before vs. After Surface Analysis**")
        comp_data = {
            "Metric": ["Surface Atoms", f"{A} on Surface", f"Surface {A} Ratio"],
            "Initial": [init_surf, init_surf_A, f"{init_ratio:.4f}"],
            "Final": [final_surf, final_surf_A, f"{final_ratio:.4f}"]
        }
        st.table(pd.DataFrame(comp_data))

        st.markdown("---")
        st.markdown("#### **Final State Sanity Checks**")
        errors = 0
        if len(res['particle']) != len(read(res['initial_xyz'])):
            st.error(f"❌ Atom count mismatch!")
            errors += 1
        else:
            st.success(f"✅ Atom count is correct: {len(res['particle'])}")
        
        final_energy = df_log['Energy (eV)'].iloc[-1] if not df_log.empty else 0
        if math.isnan(final_energy):
            st.error("❌ Final energy is NaN.")
            errors += 1
        else:
            st.success(f"✅ Final energy is valid: {final_energy:.4f} eV")

        if final_surf == 0:
            st.warning("⚠️ No surface atoms detected.")
        elif final_ratio == 0 or final_ratio == 1:
            st.warning("⚠️ Surface composition is fully segregated (0% or 100% of one element).")
        else:
            st.success("✅ Surface composition is mixed.")

    # --- 2. 3D Structures ---
    st.subheader("🔬 3D Structures")
    col1, col2 = st.columns(2)
    def plot_static_3d(xyz_path, title):
        atoms = read(xyz_path)
        fig, ax = plt.subplots(figsize=(5,5))
        plot_atoms(atoms, ax, radii=0.7, rotation=('30x,30y,0z'))
        ax.set_title(title)
        ax.axis('off')
        return fig
    
    with col1:
        st.pyplot(plot_static_3d(res['initial_xyz'], "Initial Structure"))
    with col2:
        st.pyplot(plot_static_3d(res['final_xyz'], "Final Structure"))
        
    # --- 3. Evolution Plots ---
    st.subheader("📈 Evolution Plots")
    if not df_log.empty:
        fig1, ax1 = plt.subplots()
        ax1.plot(df_log["Step"], df_log["Energy (eV)"], label="Energy")
        ax1.set(xlabel="MC Step", ylabel="Energy (eV)", title="Energy vs. Step")
        ax1.grid(True)
        st.pyplot(fig1)

        fig2, ax2 = plt.subplots()
        ax2.plot(df_log["Step"], df_log[f"Surface {A} Ratio"], color="orange")
        ax2.set(xlabel="MC Step", ylabel=f"Surface {A} Ratio", title=f"Surface {A} Ratio vs. Step")
        ax2.grid(True)
        st.pyplot(fig2)
    else:
        st.info("No log entries to plot.")

    # --- 4. Downloads & Data ---
    st.subheader("📦 Artifacts & Data")
    
    def make_download_link(path, label):
        with open(path, "rb") as f:
            b64 = base64.b64encode(f.read()).decode()
        href = f'<a href="data:application/octet-stream;base64,{b64}" download="{os.path.basename(path)}">📥 Download {label}</a>'
        return href
        
    dl_col1, dl_col2 = st.columns(2)
    with dl_col1:
        st.markdown(make_download_link(res["initial_xyz"], "Initial .xyz"), unsafe_allow_html=True)
        st.markdown(make_download_link(res["final_xyz"], "Final .xyz"), unsafe_allow_html=True)
    with dl_col2:
        st.markdown(make_download_link(res["xlsx_file"], "Log .xlsx"), unsafe_allow_html=True)
        if "video_file" in res and res["video_file"]:
            st.markdown(make_download_link(res["video_file"], "Evolution .mp4"), unsafe_allow_html=True)

    with st.expander("Show Raw Log Data"):
        st.dataframe(df_log)