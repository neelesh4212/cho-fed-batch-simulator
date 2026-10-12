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
st.caption("Educational CHO-inspired model · Illustrative parameters · Not experimentally validated")

with st.sidebar:
    st.header("Run settings")
    with st.form("settings"):
        duration = st.number_input("Culture duration (days)", 1, 30, 14)
        initial_glucose = st.number_input("Initial glucose (g/L)", 0.0, 50.0, 8.0, 1.0)
        initial_cells = st.number_input("Initial viable cells (million/mL)", 0.1, 10.0, 1.0, 0.1)
        st.subheader("Feeding")
        feed_glucose = st.number_input("Feed glucose (g/L)", 0.0, 300.0, 100.0, 10.0)
        flow = st.number_input("Feed rate (L/day)", 0.0, 0.5, 0.1, 0.05, format="%.2f")
        target = st.number_input("Glucose target (g/L)", 0.0, 30.0, 10.0, 1.0)
        volume = st.number_input("Initial volume (L)", 0.1, 2.0, 1.0, 0.1)
        cutoff = st.number_input("Feed cutoff volume (L)", 0.1, 2.0, 1.5, 0.1)
        submitted = st.form_submit_button("Run simulation", type="primary", width="stretch")

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
    ["Final product", "Final viability", "Final lactate", "Final volume"],
    ["Final_Product", "Final_Viability", "Final_Lactate", "Final_Volume"], ["g/L", "%", "g/L", "L"]):
    col.metric(label, f"{summary[key]:.2f} {unit}")
st.caption(f"Feed delivered: {summary['Total_Feed']:.2f} L · Biological step: 3 hours · Controller interval: 1 day")


def chart(title, series, ylabel, limits=None, threshold=None):
    fig = go.Figure()
    for column, label, color in series:
        fig.add_trace(go.Scatter(x=data.Days, y=data[column], name=label, mode="lines", line=dict(color=color, width=2), hovertemplate="Day %{x:.3f}<br>%{y:.3f}<extra>%{fullData.name}</extra>"))
    if threshold is not None:
        fig.add_hline(y=threshold, line_dash="dot", line_color="#777777", annotation_text="Feed cutoff")
    fig.update_layout(title=title, height=300, margin=dict(l=10,r=10,t=45,b=10), template="plotly_white", xaxis_title="Time (days)", yaxis_title=ylabel, legend=dict(orientation="h", y=1.15))
    if limits:
        fig.update_yaxes(range=limits)
    return fig


dynamics, biology, table = st.tabs(["Process dynamics", "Cell biology", "Data & downloads"])
with dynamics:
    left, right = st.columns(2)
    left.plotly_chart(chart("Product concentration", [("Product", "Product", "#008c95")], "Product (g/L)"), width="stretch")
    fig = make_subplots(specs=[[{"secondary_y": True}]])
    # Explicit interval widths distinguish a held flow from instantaneous pulses.
    fig.add_trace(go.Bar(x=(data.Days.iloc[:-1].to_numpy()+data.Days.iloc[1:].to_numpy())/2,
                        y=data.Feed_Rate.iloc[:-1], width=data.Days.diff().iloc[1:],
                        name="Feed rate", marker_color="#e99a26", opacity=0.65), secondary_y=False)
    fig.add_trace(go.Scatter(x=data.Days, y=data.Glucose, name="Glucose", mode="lines", line_color="#2878b5"), secondary_y=True)
    fig.update_layout(title="Glucose and feeding", template="plotly_white", height=300, margin=dict(l=10,r=10,t=45,b=10), legend=dict(orientation="h", y=1.15), xaxis_title="Time (days)")
    fig.update_yaxes(title_text="Feed (L/day)", rangemode="tozero", secondary_y=False)
    fig.update_yaxes(title_text="Glucose (g/L)", rangemode="tozero", secondary_y=True)
    right.plotly_chart(fig, width="stretch")
    left.plotly_chart(chart("Cell viability", [("Viability", "Viability", "#26965b")], "Viability (%)", [0,100]), width="stretch")
    right.plotly_chart(chart("Reactor volume", [("Reactor_Volume", "Volume", "#e99a26")], "Volume (L)", threshold=parameters["max_feed_volume"]), width="stretch")
with biology:
    left, right = st.columns(2)
    left.plotly_chart(chart("Cell populations", [("Viable_Cells","Viable","#26965b"),("Dead_Cells","Dead","#c65353")], "Million cells/mL"), width="stretch")
    right.plotly_chart(chart("Metabolites", [("Lactate","Lactate","#e99a26"),("Ammonia","Ammonia","#9467bd")], "Concentration (g/L)"), width="stretch")
    left.plotly_chart(chart("Growth and death", [("Growth_Rate","Growth","#008c95"),("Death_Rate","Death","#c65353")], "Rate (1/day)"), width="stretch")
    right.plotly_chart(chart("Lactate metabolism", [("Dynamic_qL","Net lactate rate","#2878b5")], "pg/(cell.day)"), width="stretch")
with table:
    st.download_button("Download results CSV", data.to_csv(index=False).encode(), "simulation.csv", "text/csv")
    st.download_button("Download parameters", json.dumps(parameters, indent=2), "parameters.json", "application/json")
    st.dataframe(data, hide_index=True, width="stretch")
with st.expander("Model assumptions"):
    st.markdown("Rates are held constant within each three-hour step. Metabolites use the starting viable-cell population. Feeding decisions are daily. No oxygen, pH, glutamine, sampling losses or dead-cell lysis are modeled. Glucose uptake is capped, but growth is not recalculated at an intra-step depletion event. Runs are held in session memory; no shared database is updated.")
