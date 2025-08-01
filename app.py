# --- Imports ---
import streamlit as st
import numpy as np
import pandas as pd
import os
import threading
import time
from ase.build import bulk
from ase.io import write, read
from ase import Atoms
from ase.neighborlist import neighbor_list
from openpyxl import Workbook
import matplotlib.pyplot as plt
import tempfile
import zipfile
from io import BytesIO

# --- Streamlit Config ---
st.set_page_config(
    page_title="Monte Carlo Nanoparticle Simulator",
    layout="wide",
    page_icon="🧪"
)
st.markdown("<h1 style='text-align: center;'>🧬 Monte Carlo Surface Segregation Simulator</h1>", unsafe_allow_html=True)
st.markdown("---")

# --- Initialize Session State ---
if 'simulation_running' not in st.session_state:
    st.session_state.simulation_running = False

if 'result' not in st.session_state:
    st.session_state.result = None

if 'progress' not in st.session_state:
    st.session_state.progress = 0

if 'status_text' not in st.session_state:
    st.session_state.status_text = "Waiting to start..."
# --- Sidebar UI: Scientifically Accurate Input Panel ---
with st.sidebar:
    st.markdown("## ⚙️ Simulation Configuration")

    # --- Element Definitions ---
    st.markdown("### 🧪 Elements & Composition")
    element_a = st.text_input("Element A Symbol", value="Pt", help="Primary component (e.g. Pt)")
    element_b = st.text_input("Element B Symbol", value="Ru", help="Secondary component (e.g. Ru)")
    composition_a = st.slider(
        f"Mole Fraction of {element_a} (%)",
        min_value=0, max_value=100, value=50,
        help="Percentage of Element A atoms in the total composition"
    )

    # --- Thermodynamic Controls ---
    st.markdown("### 🌡️ Thermodynamic Parameters")
    temperature = st.number_input(
        "Simulation Temperature (K)",
        min_value=1, max_value=5000, value=300, step=10,
        help="Thermal energy of system (Kelvin)"
    )
    mc_steps = st.number_input(
        "Monte Carlo Steps",
        min_value=1000, max_value=1_000_000, value=10000, step=1000,
        help="Total number of Monte Carlo iterations"
    )

    # --- Lattice Configuration ---
    st.markdown("### 🧱 Nanoparticle Structure")
    lattice_type = st.selectbox("Crystal Lattice", ["fcc", "bcc", "hcp"], help="Base lattice type for nanoparticle")
    surface_facet = st.selectbox("Surface Orientation", ["(111)", "(100)", "(110)", "(001)"], help="Low-index surface type")

    col1, col2, col3 = st.columns(3)
    layers_x = col1.slider("X Layers", 1, 20, 7)
    layers_y = col2.slider("Y Layers", 1, 20, 7)
    layers_z = col3.slider("Z Layers", 1, 20, 7)

    # --- Interaction Energies ---
    st.markdown("### 🧠 Interaction Energies (in eV)")
    st.caption("Interaction coefficients used for Metropolis energy calculations")

    with st.expander("🔗 Atom-Atom Interaction Energies"):
        col1, col2 = st.columns(2)
        coeff_aa = col1.number_input(f"{element_a}–{element_a} (ε<sub>AA</sub>)", value=-0.500, step=0.01, format="%.3f")
        coeff_bb = col1.number_input(f"{element_b}–{element_b} (ε<sub>BB</sub>)", value=-0.300, step=0.01, format="%.3f")
        coeff_ab = col2.number_input(f"{element_a}–{element_b} (ε<sub>AB</sub>)", value=-0.400, step=0.01, format="%.3f")

    with st.expander("📉 Surface Interaction Energies"):
        col1, col2 = st.columns(2)
        coeff_a_surf = col1.number_input(f"{element_a}–Surface (ε<sub>A–surf</sub>)", value=-0.200, step=0.01, format="%.3f")
        coeff_b_surf = col2.number_input(f"{element_b}–Surface (ε<sub>B–surf</sub>)", value=-0.100, step=0.01, format="%.3f")
        coeff_out_plane = st.number_input("Out-of-Plane Penalty (ΔE<sub>⊥</sub>)", value=0.150, step=0.01, format="%.3f")

    # --- Output Controls ---
    st.markdown("### 💾 Output Settings")
    snapshot_interval = st.number_input(
        "Save Interval (MC Steps)",
        min_value=10, max_value=100000, value=500, step=100,
        help="Save snapshot every N steps"
    )
    save_excel = st.checkbox("📊 Export Energy Log (Excel)", value=True)

    st.markdown("---")
    run_sim = st.button("🚀 Run Monte Carlo Simulation", use_container_width=True)
	# --- Phase 3: Build Config from UI + Scientific Presets ---

def load_preset(preset_name):
    """Return preset interaction coefficients for common bimetallic systems."""
    presets = {
        "Pt–Ru": {
            "element_a": "Pt", "element_b": "Ru",
            "coeff_aa": -0.53, "coeff_bb": -0.30, "coeff_ab": -0.41,
            "coeff_a_surf": -0.25, "coeff_b_surf": -0.18, "coeff_out_plane": 0.18
        },
        "Pt–Co": {
            "element_a": "Pt", "element_b": "Co",
            "coeff_aa": -0.52, "coeff_bb": -0.35, "coeff_ab": -0.43,
            "coeff_a_surf": -0.26, "coeff_b_surf": -0.20, "coeff_out_plane": 0.20
        },
        "Fe–Sn": {
            "element_a": "Fe", "element_b": "Sn",
            "coeff_aa": -0.48, "coeff_bb": -0.22, "coeff_ab": -0.31,
            "coeff_a_surf": -0.23, "coeff_b_surf": -0.17, "coeff_out_plane": 0.14
        },
    }
    return presets.get(preset_name)

# --- Optional Preset Selector ---
with st.sidebar.expander("📘 Load Scientific Presets"):
    preset_choice = st.selectbox("Preset System", ["None", "Pt–Ru", "Pt–Co", "Fe–Sn"], index=0)
    if preset_choice != "None":
        preset = load_preset(preset_choice)
        if preset:
            element_a = preset["element_a"]
            element_b = preset["element_b"]
            coeff_aa = preset["coeff_aa"]
            coeff_bb = preset["coeff_bb"]
            coeff_ab = preset["coeff_ab"]
            coeff_a_surf = preset["coeff_a_surf"]
            coeff_b_surf = preset["coeff_b_surf"]
            coeff_out_plane = preset["coeff_out_plane"]
            st.success(f"{preset_choice} parameters loaded.")

# --- Build Config Dictionary ---
simulation_config = {
    "element_a": element_a,
    "element_b": element_b,
    "composition_a": composition_a / 100.0,  # Convert % to fraction
    "temperature": temperature,
    "mc_steps": mc_steps,
    "lattice_type": lattice_type,
    "surface_facet": surface_facet,
    "layers": {
        "x": layers_x,
        "y": layers_y,
        "z": layers_z,
    },
    "coefficients": {
        f"{element_a}-{element_a}": coeff_aa,
        f"{element_b}-{element_b}": coeff_bb,
        f"{element_a}-{element_b}": coeff_ab,
        f"{element_a}-S": coeff_a_surf,
        f"{element_b}-S": coeff_b_surf,
        "out-plane": coeff_out_plane,
    },
    "snapshot_interval": snapshot_interval,
    "save_excel": save_excel
}
from ase.lattice.cubic import FaceCenteredCubic, BodyCenteredCubic
from ase.lattice.hexagonal import HexagonalClosedPacked
from ase import Atoms

def get_base_lattice(element, lattice_type):
    """Return ASE unit cell for a single-element lattice"""
    if lattice_type == "fcc":
        return FaceCenteredCubic(symbol=element, size=(1, 1, 1))
    elif lattice_type == "bcc":
        return BodyCenteredCubic(symbol=element, size=(1, 1, 1))
    elif lattice_type == "hcp":
        return HexagonalClosedPacked(symbol=element, size=(1, 1, 1))
    else:
        raise ValueError("Unsupported lattice type")

def generate_mixed_nanoparticle(cfg):
    """Generate the nanoparticle structure using ASE with A/B atoms."""
    total_atoms = cfg["layers"]["x"] * cfg["layers"]["y"] * cfg["layers"]["z"] * 4  # 4 atoms per FCC unit cell

    num_a = int(cfg["composition_a"] * total_atoms)
    num_b = total_atoms - num_a

    atoms = get_base_lattice("X", cfg["lattice_type"]) * (cfg["layers"]["x"], cfg["layers"]["y"], cfg["layers"]["z"])
    positions = atoms.get_positions()

    symbols = [cfg["element_a"]] * num_a + [cfg["element_b"]] * num_b
    np.random.shuffle(symbols)

    mixed_atoms = Atoms(symbols=symbols, positions=positions, cell=atoms.get_cell(), pbc=True)

    return mixed_atoms, num_a, num_b
import random
import math
import time
import numpy as np
from ase.neighborlist import NeighborList

def symbol_type(sym, element_a):
    """Return 'A' if symbol matches element_a, else 'B'"""
    return 'A' if sym == element_a else 'B'

def calculate_energy(atoms, coeffs, bulk_coord=12, element_a=None):
    """
    Calculate total energy of the nanoparticle configuration.

    Energy is based on pair interactions, surface atoms have fewer neighbors.

    Args:
        atoms (ase.Atoms): nanoparticle structure
        coeffs (dict): interaction coefficients like {'xA-A':..., 'xA-S':...}
        bulk_coord (int): coordination number in bulk (e.g. 12 for FCC)
        element_a (str): symbol for atom type A (e.g., 'Pt')

    Returns:
        float: total energy
    """
    # Build neighbor list with cutoff radius ~3.0 Å (typical metallic bonding)
    cutoff = 3.0
    cutoffs = [cutoff] * len(atoms)
    nl = NeighborList(cutoffs, self_interaction=False, bothways=True)
    nl.update(atoms)

    counts = {key: 0 for key in coeffs.keys()}

    for i in range(len(atoms)):
        t_i = symbol_type(atoms[i].symbol, element_a)
        indices, offsets = nl.get_neighbors(i)
        n_neighbors = len(indices)

        # Surface atom detection: neighbors < bulk_coord
        if n_neighbors < bulk_coord:
            counts[f'x{t_i}-S'] += 1

        # Count pairs: only count pairs i < j once
        for idx_j in indices:
            if i < idx_j:
                t_j = symbol_type(atoms[idx_j].symbol, element_a)
                if t_i == t_j:
                    counts[f'x{t_i}-{t_j}'] += 1
                else:
                    counts['xA-B'] += 1

    energy = 0.0
    for key, val in counts.items():
        energy += coeffs.get(key, 0.0) * val
    return energy

def monte_carlo_simulation(atoms, coeffs, element_a, temperature, n_steps, save_interval, bulk_coord=12):
    """
    Run the Monte Carlo atom swapping simulation.

    Args:
        atoms (ase.Atoms): initial structure
        coeffs (dict): energy coefficients
        element_a (str): element symbol for A
        temperature (float): temperature in K
        n_steps (int): number of MC steps
        save_interval (int): steps between data logs
        bulk_coord (int): coordination number for bulk atoms

    Returns:
        tuple: (final_atoms, log_list)
            final_atoms (ase.Atoms): structure after MC
            log_list (list of dicts): step-wise logs of energy, surface ratios etc.
    """
    BOLTZMANN_K = 8.617333262e-5
    energy = calculate_energy(atoms, coeffs, bulk_coord, element_a)
    n_atoms = len(atoms)

    log = []
    start_time = time.time()

    for step in range(1, n_steps + 1):
        i = random.randint(0, n_atoms - 1)

        # Build neighbor list for selected atom i
        cutoff = 3.0
        cutoffs = [cutoff] * n_atoms
        nl = NeighborList(cutoffs, self_interaction=False, bothways=True)
        nl.update(atoms)
        neighbors, _ = nl.get_neighbors(i)
        if len(neighbors) == 0:
            continue

        j = random.choice(neighbors)

        # Only swap if different types
        if atoms[i].symbol == atoms[j].symbol:
            continue

        # Make trial swap
        trial_atoms = atoms.copy()
        trial_atoms[i].symbol, trial_atoms[j].symbol = trial_atoms[j].symbol, trial_atoms[i].symbol

        # Calculate new energy
        new_energy = calculate_energy(trial_atoms, coeffs, bulk_coord, element_a)
        dE = new_energy - energy

        # Metropolis acceptance
        if dE < 0 or random.random() < math.exp(-dE / (BOLTZMANN_K * temperature)):
            atoms = trial_atoms
            energy = new_energy

        # Logging
        if step % save_interval == 0:
            total_a = sum(1 for a in atoms if a.symbol == element_a)
            nl.update(atoms)
            surf_atoms = 0
            surf_a = 0
            for idx in range(n_atoms):
                neighbors, _ = nl.get_neighbors(idx)
                if len(neighbors) < bulk_coord:
                    surf_atoms += 1
                    if atoms[idx].symbol == element_a:
                        surf_a += 1
            surf_ratio = surf_a / surf_atoms if surf_atoms > 0 else 0

            log.append({
                'Step': step,
                'Energy (eV)': energy,
                f'Total {element_a}': total_a,
                'Surface Atoms': surf_atoms,
                f'{element_a} on Surface': surf_a,
                f'Surface {element_a} Ratio': surf_ratio
            })

    total_time = time.time() - start_time
    print(f"Monte Carlo simulation completed in {total_time:.2f} seconds")

    return atoms, log
import pandas as pd
import matplotlib.pyplot as plt
import streamlit as st

def plot_simulation_log(log, element_a):
    """
    Plot Energy and Surface Ratio vs MC Step using matplotlib.

    Args:
        log (list of dicts): simulation log returned by monte_carlo_simulation()
        element_a (str): element symbol for atom type A
    """
    if not log:
        st.warning("No simulation data to plot yet.")
        return

    # Convert log to DataFrame for easier plotting
    df = pd.DataFrame(log)

    fig, axs = plt.subplots(1, 2, figsize=(12, 5), dpi=100)

    # Energy vs Step
    axs[0].plot(df['Step'], df['Energy (eV)'], color='blue', marker='o', markersize=4, linestyle='-')
    axs[0].set_xlabel('MC Step')
    axs[0].set_ylabel('Energy (eV)')
    axs[0].set_title('Energy vs Monte Carlo Step')
    axs[0].grid(True)

    # Surface Ratio vs Step
    surface_ratio_key = f'Surface {element_a} Ratio'
    if surface_ratio_key in df.columns:
        axs[1].plot(df['Step'], df[surface_ratio_key], color='green', marker='o', markersize=4, linestyle='-')
        axs[1].set_xlabel('MC Step')
        axs[1].set_ylabel(f'Surface {element_a} Ratio')
        axs[1].set_title(f'Surface {element_a} Ratio vs MC Step')
        axs[1].grid(True)

    plt.tight_layout()
    st.pyplot(fig)

def streamlit_simulation_visualization(log, element_a):
    """
    Streamlit section to display simulation logs and plots interactively.
    """
    st.subheader("🔬 Simulation Data Logs")
    if log:
        df_log = pd.DataFrame(log)
        st.dataframe(df_log)
    else:
        st.info("Simulation log is empty. Run simulation first.")

    st.subheader("📊 Simulation Metrics Plots")
    plot_simulation_log(log, element_a)
import io
import pandas as pd
from ase.io import write
import streamlit as st

def export_simulation_log_to_excel(log, simulation_params, filename="simulation_log.xlsx"):
    """
    Export the simulation log and parameters to an Excel file in-memory.
    
    Args:
        log (list of dicts): simulation log data
        simulation_params (dict): simulation input parameters (for metadata sheet)
        filename (str): default filename for the download button

    Returns:
        bytes: Excel file in bytes buffer
    """
    output = io.BytesIO()
    with pd.ExcelWriter(output, engine='openpyxl') as writer:
        # Write simulation params in first sheet
        params_df = pd.DataFrame(list(simulation_params.items()), columns=['Parameter', 'Value'])
        params_df.to_excel(writer, sheet_name='Parameters', index=False)

        # Write log data in second sheet
        log_df = pd.DataFrame(log)
        log_df.to_excel(writer, sheet_name='Simulation Log', index=False)
    
    return output.getvalue()

def export_final_structure_to_xyz(atoms, filename="final_structure.xyz"):
    """
    Export ASE Atoms object to an XYZ file bytes buffer.
    
    Args:
        atoms (ase.Atoms): final nanoparticle structure
    
    Returns:
        bytes: XYZ file in bytes buffer
    """
    output = io.StringIO()
    write(output, atoms, format='xyz')
    xyz_data = output.getvalue()
    return xyz_data.encode('utf-8')

def streamlit_export_section(log, params, atoms):
    """
    Streamlit UI section to export and download log & structure files.
    """
    st.subheader("💾 Export & Download Results")

    # Export Excel log file
    excel_bytes = export_simulation_log_to_excel(log, params)
    st.download_button(
        label="Download Simulation Log (Excel)",
        data=excel_bytes,
        file_name="simulation_log.xlsx",
        mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
    )

    # Export XYZ file for final structure
    xyz_bytes = export_final_structure_to_xyz(atoms)
    st.download_button(
        label="Download Final Nanoparticle Structure (XYZ)",
        data=xyz_bytes,
        file_name="final_structure.xyz",
        mime="chemical/x-xyz"
    )
import py3Dmol
import streamlit as st
import streamlit.components.v1 as components

def view_structure_3d(xyz_str: str, atom_colors=None, sphere_radius=1.5, width=700, height=500):
    """
    Generates a 3Dmol.js viewer HTML string for an XYZ structure.

    Args:
        xyz_str (str): contents of an XYZ file (string)
        atom_colors (dict): mapping atom symbol -> color string, e.g. {'Pt':'gold', 'Ru':'green'}
        sphere_radius (float): radius of atom spheres
        width (int): viewer width in pixels
        height (int): viewer height in pixels

    Returns:
        None (renders in Streamlit)
    """
    view = py3Dmol.view(width=width, height=height)
    view.addModel(xyz_str, 'xyz')

    # Default colors if none provided
    if atom_colors is None:
        atom_colors = {}

    # Apply styles per atom symbol
    for atom_symbol, color in atom_colors.items():
        view.setStyle({'elem': atom_symbol}, {'sphere': {'color': color, 'radius': sphere_radius}})

    view.setBackgroundColor('white')
    view.zoomTo()
    html_str = view._make_html()

    # Render HTML with Streamlit components
    components.html(html_str, height=height + 30)
st.subheader("🔬 3D Visualization of Nanoparticle Structure")

# Load or generate initial and final XYZ strings
with open(initial_xyz_file, 'r') as f:
    initial_xyz_str = f.read()

with open(final_xyz_file, 'r') as f:
    final_xyz_str = f.read()

atom_colors = {
    user_input['element_A']: 'gold',
    user_input['element_B']: 'green',
}

st.markdown("**Initial Structure:**")
view_structure_3d(initial_xyz_str, atom_colors=atom_colors)

st.markdown("**Final Structure:**")
view_structure_3d(final_xyz_str, atom_colors=atom_colors)
import plotly.graph_objs as go
import streamlit as st

# Initialize empty lists to store simulation data
steps_data = []
energy_data = []
surface_ratio_data = []

# Create Streamlit placeholders for plots
energy_plot_placeholder = st.empty()
surface_plot_placeholder = st.empty()
def update_plots(steps, energy, surface_ratio, element_symbol):
    # Energy plot
    fig_energy = go.Figure()
    fig_energy.add_trace(go.Scatter(x=steps, y=energy, mode='lines+markers', name='Energy (eV)', line=dict(color='blue')))
    fig_energy.update_layout(title="Energy vs Monte Carlo Step",
                             xaxis_title="MC Step",
                             yaxis_title="Energy (eV)",
                             template="plotly_white")

    # Surface composition plot
    fig_surface = go.Figure()
    fig_surface.add_trace(go.Scatter(x=steps, y=surface_ratio, mode='lines+markers',
                                     name=f'Surface {element_symbol} Ratio',
                                     line=dict(color='green')))
    fig_surface.update_layout(title=f"Surface {element_symbol} Composition Ratio vs MC Step",
                              xaxis_title="MC Step",
                              yaxis_title=f"Surface {element_symbol} Ratio",
                              template="plotly_white")

    return fig_energy, fig_surface
import time

# Inside simulation loop (pseudo-code)
for step in range(1, N_STEPS + 1):
    # ... your MC step logic ...

    if step % SAVE_INTERVAL == 0:
        # Collect data
        steps_data.append(step)
        energy_data.append(current_energy)
        surface_ratio_data.append(current_surface_ratio)

        # Update plots dynamically
        fig_energy, fig_surface = update_plots(steps_data, energy_data, surface_ratio_data, user_input['element_A'])
        energy_plot_placeholder.plotly_chart(fig_energy, use_container_width=True)
        surface_plot_placeholder.plotly_chart(fig_surface, use_container_width=True)

        # Optional: small delay to allow UI refresh (adjust as needed)
        time.sleep(0.01)
import threading
import streamlit as st

if 'simulation_running' not in st.session_state:
    st.session_state.simulation_running = False
if 'simulation_thread' not in st.session_state:
    st.session_state.simulation_thread = None
if 'simulation_data' not in st.session_state:
    st.session_state.simulation_data = {
        'steps': [],
        'energies': [],
        'surface_ratios': []
    }
def run_simulation(user_input):
    # Reset data at start
    st.session_state.simulation_data = {
        'steps': [],
        'energies': [],
        'surface_ratios': []
    }
    N_STEPS = user_input['n_steps']
    
    # Initialize particle and energy here...

    for step in range(1, N_STEPS + 1):
        if not st.session_state.simulation_running:
            break  # Exit if stopped

        # Monte Carlo simulation step logic here
        # ...

        if step % user_input['save_interval'] == 0:
            # Append data safely
            st.session_state.simulation_data['steps'].append(step)
            st.session_state.simulation_data['energies'].append(current_energy)
            st.session_state.simulation_data['surface_ratios'].append(current_surface_ratio)
        
        # Allow UI update
        time.sleep(0.01)
    
    st.session_state.simulation_running = False
start_button = st.button("Start Simulation")
stop_button = st.button("Stop Simulation")

if start_button and not st.session_state.simulation_running:
    st.session_state.simulation_running = True
    sim_thread = threading.Thread(target=run_simulation, args=(user_input,), daemon=True)
    st.session_state.simulation_thread = sim_thread
    sim_thread.start()

if stop_button and st.session_state.simulation_running:
    st.session_state.simulation_running = False
    st.write("Simulation stopping... please wait.")
import plotly.graph_objs as go

energy_plot = st.empty()
surface_plot = st.empty()

def plot_data():
    steps = st.session_state.simulation_data['steps']
    energies = st.session_state.simulation_data['energies']
    surface_ratios = st.session_state.simulation_data['surface_ratios']
    element_symbol = user_input['element_A']

    fig_energy = go.Figure(data=[go.Scatter(x=steps, y=energies, mode='lines+markers')])
    fig_energy.update_layout(title='Energy vs MC Step', xaxis_title='Step', yaxis_title='Energy (eV)', template='plotly_white')

    fig_surface = go.Figure(data=[go.Scatter(x=steps, y=surface_ratios, mode='lines+markers')])
    fig_surface.update_layout(title=f'Surface {element_symbol} Ratio vs MC Step', xaxis_title='Step', yaxis_title=f'Surface {element_symbol} Ratio', template='plotly_white')

    energy_plot.plotly_chart(fig_energy, use_container_width=True)
    surface_plot.plotly_chart(fig_surface, use_container_width=True)

# Poll data every 500 ms
import time
while st.session_state.simulation_running:
    plot_data()
    time.sleep(0.5)

# Plot final data after simulation ends
plot_data()
import io

# After simulation finishes, final_xyz is path to final structure file

with open(final_xyz, 'r') as f:
    xyz_data = f.read()

st.download_button(
    label="Download Final Nanoparticle Structure (.xyz)",
    data=xyz_data,
    file_name=final_xyz,
    mime="chemical/x-xyz"
)
import openpyxl

# Save Excel log in-memory for download
import io

xlsx_path = "MMC_log.xlsx"  # your saved Excel file path

with open(xlsx_path, "rb") as f:
    xlsx_bytes = f.read()

st.download_button(
    label="Download Simulation Log (.xlsx)",
    data=xlsx_bytes,
    file_name=xlsx_path,
    mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
)
import pandas as pd

# Suppose log is a list of dicts (simulation log)
df_log = pd.DataFrame(log)  # your simulation log data

csv_buffer = io.StringIO()
df_log.to_csv(csv_buffer, index=False)
csv_data = csv_buffer.getvalue()

st.download_button(
    label="Download Simulation Log (.csv)",
    data=csv_data,
    file_name="MMC_log.csv",
    mime="text/csv"
)
import zipfile
import os

zip_path = "snapshots.zip"
snapshot_folder = "trajectory"

with zipfile.ZipFile(zip_path, 'w') as z:
    for filename in os.listdir(snapshot_folder):
        if filename.endswith(".xyz"):
            z.write(os.path.join(snapshot_folder, filename), arcname=filename)

with open(zip_path, 'rb') as f:
    zip_bytes = f.read()

st.download_button(
    label="Download All Snapshots (.zip)",
    data=zip_bytes,
    file_name=zip_path,
    mime="application/zip"
)
import io
import matplotlib.pyplot as plt

fig, ax = plt.subplots()
ax.plot([1, 2, 3], [4, 5, 6])
ax.set_title("Example Plot")

buf = io.BytesIO()
fig.savefig(buf, format='png')
buf.seek(0)

st.download_button(
    label="Download Plot Image (.png)",
    data=buf,
    file_name="plot.png",
    mime="image/png"
)
video_path = "simulation_movie.mp4"
with open(video_path, "rb") as f:
    video_bytes = f.read()

st.download_button(
    label="Download Simulation Video (.mp4)",
    data=video_bytes,
    file_name=video_path,
    mime="video/mp4"
)
import py3Dmol
import streamlit as st

def show_py3dmol_xyz(xyz_string, atom_colors):
    view = py3Dmol.view(width=400, height=400)
    view.addModel(xyz_string, 'xyz')
    for elem, color in atom_colors.items():
        view.setStyle({'elem': elem}, {'sphere': {'color': color, 'radius': 1.5}})
    view.setBackgroundColor('white')
    view.zoomTo()
    html_str = view._make_html()
    st.components.v1.html(html_str, height=450)

# Usage:
with open(final_xyz, 'r') as f:
    xyz_data = f.read()

atom_colors = {'Pt': 'silver', 'Ru': 'royalblue'}

st.markdown("### Final Nanoparticle Structure")
show_py3dmol_xyz(xyz_data, atom_colors)
video_file = open("simulation_movie.mp4", "rb").read()
st.video(video_file)
import plotly.graph_objects as go

def plot_interactive_simulation_log(log):
    steps = [entry['Step'] for entry in log]
    energy = [entry['Energy (eV)'] for entry in log]
    surface_ratio = [entry['Surface Pt Ratio'] for entry in log]

    fig = go.Figure()
    fig.add_trace(go.Scatter(x=steps, y=energy, mode='lines+markers', name='Energy (eV)', line=dict(color='blue')))
    fig.add_trace(go.Scatter(x=steps, y=surface_ratio, mode='lines+markers', name='Surface Pt Ratio', line=dict(color='green'), yaxis="y2"))

    fig.update_layout(
        title="Simulation Log",
        xaxis_title="Monte Carlo Step",
        yaxis=dict(title='Energy (eV)', side='left'),
        yaxis2=dict(title='Surface Pt Ratio', overlaying='y', side='right'),
        legend=dict(x=0.1, y=1.1, orientation='h')
    )
    st.plotly_chart(fig, use_container_width=True)

# Usage
plot_interactive_simulation_log(log)
import streamlit as st
import tempfile
import os

def export_files(particle, log, trajectory_folder):
    """
    Prepare downloadable files: final structure XYZ, Excel log, and optionally zipped snapshots.
    """

    # 1. Export final structure as XYZ file
    final_xyz_path = tempfile.NamedTemporaryFile(delete=False, suffix=".xyz").name
    write(final_xyz_path, particle)

    # 2. Export log as Excel file
    from openpyxl import Workbook

    wb = Workbook()
    ws = wb.active
    ws.title = "Simulation Log"

    # Write headers and data
    headers = list(log[0].keys()) if log else []
    for col_num, header in enumerate(headers, 1):
        ws.cell(row=1, column=col_num, value=header)

    for row_num, entry in enumerate(log, 2):
        for col_num, header in enumerate(headers, 1):
            ws.cell(row=row_num, column=col_num, value=entry[header])

    excel_path = tempfile.NamedTemporaryFile(delete=False, suffix=".xlsx").name
    wb.save(excel_path)

    # 3. Optionally create a zip of trajectory snapshots if folder exists
    import shutil
    zip_path = None
    if os.path.exists(trajectory_folder) and os.listdir(trajectory_folder):
        zip_path = tempfile.NamedTemporaryFile(delete=False, suffix=".zip").name
        shutil.make_archive(base_name=zip_path.replace('.zip', ''),
                            format='zip',
                            root_dir=trajectory_folder)

    return final_xyz_path, excel_path, zip_path

def download_section(particle, log, trajectory_folder="trajectory"):
    st.subheader("📁 Download Simulation Outputs")

    xyz_path, excel_path, zip_path = export_files(particle, log, trajectory_folder)

    st.markdown("**Final Nanoparticle Structure (.xyz):**")
    with open(xyz_path, "rb") as f:
        st.download_button(
            label="Download XYZ",
            data=f,
            file_name=os.path.basename(xyz_path),
            mime="chemical/x-xyz"
        )

    st.markdown("**Simulation Log (.xlsx):**")
    with open(excel_path, "rb") as f:
        st.download_button(
            label="Download Excel Log",
            data=f,
            file_name=os.path.basename(excel_path),
            mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
        )

    if zip_path:
        st.markdown("**Trajectory Snapshots (.zip):**")
        with open(zip_path, "rb") as f:
            st.download_button(
                label="Download Trajectory Zip",
                data=f,
                file_name=os.path.basename(zip_path),
                mime="application/zip"
            )

    # Cleanup temporary files on app exit (optional)
    import atexit
    def cleanup_files():
        for p in [xyz_path, excel_path, zip_path]:
            if p and os.path.exists(p):
                os.remove(p)
    atexit.register(cleanup_files)
import streamlit as st

def validate_inputs(user_input):
    errors = []
    if user_input['composition_A'] < 0 or user_input['composition_A'] > 1:
        errors.append("Composition of element A must be between 0 and 1.")

    if user_input['temperature'] <= 0:
        errors.append("Temperature must be greater than zero Kelvin.")

    if user_input['n_steps'] <= 0:
        errors.append("Number of Monte Carlo steps must be positive.")

    if user_input['lattice_type'].lower() not in ['fcc', 'bcc', 'hcp']:
        errors.append("Lattice type must be one of: fcc, bcc, hcp.")

    # Add more validations as needed

    return errors

def error_handling_section(user_input):
    st.subheader("⚠️ Input Validation")

    errors = validate_inputs(user_input)
    if errors:
        for err in errors:
            st.error(err)
        st.stop()
    else:
        st.success("All inputs are valid. Ready to run simulation.")
import streamlit as st
import threading

def run_simulation_in_thread(simulation_func, *args, **kwargs):
    result = {}
    def target():
        result['output'] = simulation_func(*args, **kwargs)
    thread = threading.Thread(target=target)
    thread.start()
    return thread, result

def simulation_progress_ui():
    st.info("Simulation running... Please wait.")
    # Optionally add a spinner or progress bar updating with real-time status

# Usage inside app.py:
# thread, result = run_simulation_in_thread(monte_carlo_simulation, particle, coeffs, ... )
# while thread.is_alive():
#     simulation_progress_ui()
#     time.sleep(1)
# # After completion:
# particle, log = result['output']

