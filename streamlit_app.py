import json
import sys
from pathlib import Path
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import streamlit as st

sys.path.insert(0, str(Path(__file__).resolve().parent / "src"))
from web_simulation import simulate

st.set_page_config(page_title="CHO Fed-Batch Simulator", page_icon=":test_tube:", layout="wide")
st.title("CHO Fed-Batch Simulator")
st.caption("Educational CHO-Inspired Model · Illustrative Parameters")
st.markdown("""
<style>
[data-testid="stMetricLabel"] p { font-size: 18px; font-weight: 600; }
[data-testid="stTabs"] [role="tab"] p { font-size: 18px; }
[data-testid="stTabs"] [role="tablist"] { flex-wrap: wrap; height: auto; }
</style>
""", unsafe_allow_html=True)

with st.sidebar:
    st.header("Reactor Settings")
    with st.form("settings"):
        duration = st.number_input("Culture Duration (days)", 1, 30, 14)
        initial_glucose = st.number_input("Initial Glucose (g/L)", 0.0, 50.0, 8.0, 1.0)
        initial_cells = st.number_input("Initial Viable Cells (million/mL)", 0.1, 10.0, 1.0, 0.1)
        st.subheader("Feeding")
        feed_glucose = st.number_input("Feed Glucose (g/L)", 0.0, 300.0, 100.0, 10.0)
        flow = st.number_input("Feed Rate (L/day)", 0.0, 0.5, 0.1, 0.05, format="%.2f")
        target = st.number_input("Glucose Target (g/L)", 0.0, 30.0, 10.0, 1.0)
        volume = st.number_input("Initial Volume (L)", 0.1, 2.0, 1.0, 0.1)
        cutoff = st.number_input("Feed Cutoff Volume (L)", 0.1, 2.0, 1.5, 0.1)
        submitted = st.form_submit_button("Run Simulation", type="primary", width="stretch")

inputs = dict(simulation_duration=float(duration), initial_glucose=initial_glucose,
              initial_cells=initial_cells, feed_glucose=feed_glucose, feed_rate_on=flow,
              glucose_target=target, initial_volume=volume, max_feed_volume=cutoff)
if submitted or "run" not in st.session_state:
    try:
        st.session_state.run = simulate(inputs)
    except ValueError as error:
        st.error(str(error))
        st.stop()

data, summary, parameters = st.session_state.run
if inputs != {key: parameters[key] for key in inputs}:
    st.info("Settings differ from the displayed run.")
metrics = st.columns(4)
for col, label, key, unit in zip(metrics,
    ["Final Product", "Final Viability", "Final Lactate", "Final Volume"],
    ["Final_Product", "Final_Viability", "Final_Lactate", "Final_Volume"], ["g/L", "%", "g/L", "L"]):
    col.metric(label, f"{summary[key]:.2f} {unit}")
st.caption(f"Feed Delivered: {summary['Total_Feed']:.2f} L · Calculation Interval: 3 hours · Controller Interval: 1 day")


def chart(title, series, ylabel, limits=None, threshold=None):
    fig = go.Figure()
    for column, label, color in series:
        fig.add_trace(go.Scatter(x=data.Days, y=data[column], name=label, mode="lines", line=dict(color=color, width=2), hovertemplate="Day %{x:.3f}<br>%{y:.3f}<extra>%{fullData.name}</extra>"))
    if threshold is not None:
        fig.add_hline(y=threshold, line_dash="dot", line_color="#777777", annotation_text="Feed Cutoff")
    # Reserve a separate header row for the legend below the chart title.
    fig.update_layout(title=dict(text=title, y=0.98, yanchor="top", font=dict(size=16)), height=330, margin=dict(l=10,r=10,t=95,b=10), template="plotly_white", xaxis_title="Time (days)", yaxis_title=ylabel, legend=dict(orientation="h", y=1.02, yanchor="bottom", x=0))
    if limits:
        fig.update_yaxes(range=limits)
    return fig


dynamics, biology, table = st.tabs(["Process Dynamics", "Cell Biology", "Data and Downloads"])
with dynamics:
    left, right = st.columns(2)
    left.plotly_chart(chart("Product Concentration", [("Product", "Product", "#008c95")], "Product (g/L)"), width="stretch")
    fig = make_subplots(specs=[[{"secondary_y": True}]])
    # Explicit interval widths distinguish a held flow from instantaneous pulses.
    fig.add_trace(go.Bar(x=(data.Days.iloc[:-1].to_numpy()+data.Days.iloc[1:].to_numpy())/2,
                        y=data.Feed_Rate.iloc[:-1], width=data.Days.diff().iloc[1:],
                        name="Feed Rate", marker_color="#e99a26", marker_line_width=0, opacity=1.0), secondary_y=False)
    fig.add_trace(go.Scatter(x=data.Days, y=data.Glucose, name="Glucose", mode="lines", line_color="#2878b5"), secondary_y=True)
    fig.update_layout(title=dict(text="Glucose and Feeding", y=0.98, yanchor="top", font=dict(size=16)), template="plotly_white", bargap=0, bargroupgap=0, height=330, margin=dict(l=10,r=10,t=95,b=10), legend=dict(orientation="h", y=1.02, yanchor="bottom", x=0), xaxis_title="Time (days)")
    fig.update_yaxes(title_text="Feed (L/day)", rangemode="tozero", secondary_y=False)
    fig.update_yaxes(title_text="Glucose (g/L)", rangemode="tozero", secondary_y=True)
    right.plotly_chart(fig, width="stretch")
    left.plotly_chart(chart("Cell Viability", [("Viability", "Viability", "#26965b")], "Viability (%)", [0,100]), width="stretch")
    right.plotly_chart(chart("Reactor Volume", [("Reactor_Volume", "Volume", "#e99a26")], "Volume (L)", threshold=parameters["max_feed_volume"]), width="stretch")
with biology:
    left, right = st.columns(2)
    left.plotly_chart(chart("Cell Populations", [("Viable_Cells","Viable","#26965b"),("Dead_Cells","Dead","#c65353")], "Million cells/mL"), width="stretch")
    right.plotly_chart(chart("Metabolites", [("Lactate","Lactate","#e99a26"),("Ammonia","Ammonia","#9467bd")], "Concentration (g/L)"), width="stretch")
    left.plotly_chart(chart("Growth and Death", [("Growth_Rate","Growth","#008c95"),("Death_Rate","Death","#c65353")], "Rate (1/day)"), width="stretch")
    right.plotly_chart(chart("Lactate Metabolism", [("Dynamic_qL","Net Lactate Rate","#2878b5")], "pg/(cell.day)"), width="stretch")
with table:
    st.download_button("Download Results CSV", data.to_csv(index=False).encode(), "simulation.csv", "text/csv")
    st.download_button("Download Parameters", json.dumps(parameters, indent=2), "parameters.json", "application/json")
    st.dataframe(data, hide_index=True, width="stretch")
with st.expander("Model Assumptions"):
    st.markdown("Parameters are illustrative and have not been experimentally validated. Rates are held constant within each three-hour step. Metabolites use the starting viable-cell population. Feeding decisions are daily. No oxygen, pH, glutamine, sampling losses or dead-cell lysis are modeled. Glucose uptake is capped, but growth is not recalculated at an intra-step depletion event. Runs are held in session memory; no shared database is updated.")
