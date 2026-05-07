import pandas as pd
import os
import re
import streamlit as st
import altair as alt
from io import StringIO

pd.set_option('display.max_columns', 20)
pd.set_option('display.max_rows', 10)

dossier = r"C:\Users\AL274877\Desktop\PhreeqC-Comsol\ex11PhreeqCv.0.1\orchestra\SNIA\outputHistory"
dossier = r"C:\Users\AL274877\Desktop\PhreeqC-Comsol\CellulesSeb\Cell2\outputSpeciation"

st.set_page_config(layout="wide")
alt.data_transformers.disable_max_rows()

def extraire_numero(nom_fichier):
    match = re.search(r'(\d+)', nom_fichier)
    return int(match.group(1)) if match else float("inf")

@st.cache_data(show_spinner="Lecture des fichiers...")
def charger_sorties(dossier):
    if not os.path.exists(dossier):
        st.error("Dossier introuvable")
        return []

    fichiers_txt = sorted(
        [f for f in os.listdir(dossier) if f.endswith(".txt")],
        key=extraire_numero
    )

    dfs = []

    for i, nom_fichier in enumerate(fichiers_txt):
        chemin = os.path.join(dossier, nom_fichier)

        df = pd.read_csv(
            chemin,
            sep=r"\s+",
            comment="%",
            header=0
        )

        if df.empty:
            continue

        df.insert(0, "t", i)
        dfs.append(df)

    return dfs


df_sortiePhreeqC = charger_sorties(dossier)

if not df_sortiePhreeqC:
    st.stop()

st.title("Profil spatial PICCTS")

index_df = st.slider(
    "Choisir le pas de temps",
    0,
    len(df_sortiePhreeqC) - 1,
    0
)

df_selected = df_sortiePhreeqC[index_df]

all_columns = [col for col in df_selected.columns if col not in ['t', 'x']]

cols_to_plot = st.multiselect(
    "Choisir les colonnes à tracer :",
    all_columns,
    default=all_columns[:1]
)

if not cols_to_plot:
    st.warning("Sélectionne au moins une espèce.")
    st.stop()

plot_data = []

use_x = "x" in df_selected.columns

for col in cols_to_plot:
    for i in range(len(df_selected)):
        plot_data.append({
            "x": df_selected["x"].iloc[i] if use_x else i,
            "value": df_selected[col].iloc[i],
            "column": col
        })

plot_df = pd.DataFrame(plot_data)

chart = alt.Chart(plot_df).mark_line().encode(
    x=alt.X("x", title="Position"),
    y=alt.Y("value", title="Valeur"),
    color="column",
    tooltip=["column", "x", "value"]
).interactive()

st.altair_chart(chart, use_container_width=True)

csv_buffer = StringIO()
plot_df.to_csv(csv_buffer, index=False)

st.download_button(
    label="Télécharger les données en CSV",
    data=csv_buffer.getvalue(),
    file_name="profil_spatial.csv",
    mime="text/csv"
)