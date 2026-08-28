import streamlit as st
import pandas as pd
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import ast
import warnings
warnings.filterwarnings('ignore')

try:
    import community as community_louvain
    LOUVAIN_OK = True
except:
    LOUVAIN_OK = False

# ═══════════════════════════════════════════════════════════
# CONFIGURATION
# ═══════════════════════════════════════════════════════════
st.set_page_config(
    page_title="Analyse Communautés PPI — GA-PPI-Net",
    page_icon="🧬",
    layout="wide"
)

st.markdown("""
<style>
    .main-title { font-size:2.2rem; font-weight:700; color:#1f4e79; text-align:center; padding:1rem 0 0.2rem 0; }
    .subtitle   { font-size:1rem; color:#555; text-align:center; margin-bottom:1.5rem; }
    .section-header { font-size:1.1rem; font-weight:700; color:#1f4e79; border-bottom:2px solid #1f4e79; padding-bottom:0.3rem; margin:1rem 0 0.7rem 0; }
    .stButton>button { background:#1f4e79; color:white; border-radius:8px; border:none; font-weight:600; }
</style>
""", unsafe_allow_html=True)

# ═══════════════════════════════════════════════════════════
# FONCTIONS
# ═══════════════════════════════════════════════════════════
@st.cache_data(show_spinner="⏳ Chargement du réseau PPI...")
def charger_reseau(fichier):
    df = pd.read_csv(fichier, sep=';')
    df.columns = ['Gene1', 'Gene2', 'Interaction', 'Similarite']
    df['Interaction'] = pd.to_numeric(df['Interaction'], errors='coerce')
    df['Similarite']  = pd.to_numeric(df['Similarite'],  errors='coerce')
    df = df.dropna()
    dct = {}
    for _, row in df.iterrows():
        dct[(row['Gene1'], row['Gene2'])] = (float(row['Interaction']), float(row['Similarite']))
    return df, dct

def charger_communautes(fichier):
    comms = []
    contenu = fichier.read().decode('utf-8')
    for ligne in contenu.strip().split('\n'):
        if ligne.strip():
            try:
                comms.append(ast.literal_eval(ligne.strip()))
            except:
                pass
    return comms

def chercher_paire(g1, g2, dct):
    if (g1, g2) in dct: return dct[(g1, g2)]
    if (g2, g1) in dct: return dct[(g2, g1)]
    return None

def construire_graphe(genes, dct, seuil=150):
    from itertools import combinations
    G = nx.Graph()
    G.add_nodes_from(genes)
    for g1, g2 in combinations(genes, 2):
        p = chercher_paire(g1, g2, dct)
        if p and p[0] >= seuil:
            G.add_edge(g1, g2, weight=p[0], interaction=p[0], similarite=p[1])
    return G

def calculer_moyennes(genes, dct):
    from itertools import combinations
    sims, ints = [], []
    total = 0
    for g1, g2 in combinations(genes, 2):
        total += 1
        p = chercher_paire(g1, g2, dct)
        if p:
            ints.append(p[0]); sims.append(p[1])
    return {
        'AVGSIM': round(np.mean(sims), 4) if sims else 0,
        'AVGInteraction': round(np.mean(ints), 2) if ints else 0,
        'paires_totales': total,
        'paires_trouvees': len(sims),
        'liste_sim': sims,
        'liste_int': ints
    }

def calculer_metriques(G):
    m = {}
    n, e = G.number_of_nodes(), G.number_of_edges()
    m['nb_noeuds'] = n; m['nb_aretes'] = e
    if n == 0: return m
    degres = dict(G.degree())
    m['degre_moyen']     = round(np.mean(list(degres.values())), 3)
    m['degre_max']       = max(degres.values())
    m['degre_min']       = min(degres.values())
    m['hub_principal']   = max(degres, key=degres.get)
    m['distribution_degres'] = degres
    m['densite']         = round(nx.density(G), 4)
    m['max_aretes_possibles'] = int(n*(n-1)/2)
    m['clustering_moyen'] = round(nx.average_clustering(G), 4)
    if nx.is_connected(G):
        m['est_connexe'] = True
        m['distance_geodesique_moy'] = round(nx.average_shortest_path_length(G), 3)
        m['diametre'] = nx.diameter(G)
    else:
        m['est_connexe'] = False
        comps = list(nx.connected_components(G))
        m['nb_composantes'] = len(comps)
        grande = G.subgraph(max(comps, key=len)).copy()
        m['distance_geodesique_moy'] = round(nx.average_shortest_path_length(grande), 3) if grande.number_of_nodes() > 1 else 0
        m['diametre'] = nx.diameter(grande) if grande.number_of_nodes() > 1 else 0
    if LOUVAIN_OK and e > 0:
        try:
            partition = community_louvain.best_partition(G)
            m['modularite'] = round(community_louvain.modularity(partition, G), 4)
            m['nb_sous_groupes'] = len(set(partition.values()))
            m['partition'] = partition
        except:
            m['modularite'] = 0; m['nb_sous_groupes'] = 1; m['partition'] = None
    else:
        m['modularite'] = 'N/A'; m['nb_sous_groupes'] = 'N/A'; m['partition'] = None
    m['score_pizzuti'] = round(e * m['degre_moyen'], 3)
    return m

def dessiner_graphe(G, met, titre):
    if G.number_of_nodes() == 0:
        st.warning("⚠️ Aucune interaction trouvée."); return
    fig, ax = plt.subplots(figsize=(12, 8))
    pos = nx.spring_layout(G, seed=42, k=2.5) if G.number_of_nodes() <= 30 else nx.kamada_kawai_layout(G)
    degres = dict(G.degree())
    max_deg = max(degres.values()) if degres else 1
    couleurs = [plt.cm.YlOrRd(degres[n] / max_deg) for n in G.nodes()]
    tailles  = [300 + degres[n] * 100 for n in G.nodes()]
    poids    = [G[u][v].get('interaction', 150) / 1000 * 3 for u, v in G.edges()]
    nx.draw_networkx_edges(G, pos, width=poids, alpha=0.35, edge_color='#4a7db5', ax=ax)
    nx.draw_networkx_nodes(G, pos, node_color=couleurs, node_size=tailles, alpha=0.9, ax=ax)
    if G.number_of_nodes() <= 50:
        nx.draw_networkx_labels(G, pos, font_size=7, font_weight='bold', ax=ax)
    hub = met.get('hub_principal', '')
    if hub and hub in pos:
        idx = list(G.nodes()).index(hub)
        nx.draw_networkx_nodes(G, pos, nodelist=[hub], node_color='red', node_size=tailles[idx]+200, alpha=1.0, ax=ax)
    p1 = mpatches.Patch(color=plt.cm.YlOrRd(0.15), label='Faible degré')
    p2 = mpatches.Patch(color='red',                label=f'Hub : {hub}')
    p3 = mpatches.Patch(color=plt.cm.YlOrRd(0.85), label='Degré élevé')
    ax.legend(handles=[p1, p3, p2], loc='upper left', fontsize=8)
    ax.set_title(f"{titre}\n{G.number_of_nodes()} gènes · {G.number_of_edges()} interactions · Hub : {hub}",
                 fontsize=12, fontweight='bold', pad=12)
    ax.axis('off')
    fig.patch.set_facecolor('#f8f9fa')
    st.pyplot(fig); plt.close()

def dessiner_histogrammes(moy):
    fig, axes = plt.subplots(1, 3, figsize=(14, 4))
    fig.suptitle('Distributions des valeurs de la communauté', fontsize=12, fontweight='bold')
    degres_vals = sorted(moy.get('degres_vals', []), reverse=True)
    if degres_vals:
        axes[0].bar(range(len(degres_vals)), degres_vals, color='#1f4e79', alpha=0.7)
        axes[0].set_title('① Distribution des degrés', fontweight='bold')
        axes[0].set_xlabel('Gènes'); axes[0].set_ylabel('Degré')
    if moy['liste_sim']:
        axes[1].hist(moy['liste_sim'], bins=20, color='#2e7d32', alpha=0.75, edgecolor='white')
        axes[1].axvline(x=moy['AVGSIM'], color='red', linestyle='--', linewidth=2, label=f"AVGSIM={moy['AVGSIM']:.3f}")
        axes[1].axvline(x=0.5, color='orange', linestyle=':', linewidth=1.5, label='Seuil ∇S=0.5')
        axes[1].set_title('Similarité GS2', fontweight='bold'); axes[1].legend(fontsize=8)
    if moy['liste_int']:
        axes[2].hist(moy['liste_int'], bins=20, color='#e65100', alpha=0.75, edgecolor='white')
        axes[2].axvline(x=moy['AVGInteraction'], color='red', linestyle='--', linewidth=2, label=f"AVGInt={moy['AVGInteraction']:.0f}")
        axes[2].set_title('Score Interaction STRING', fontweight='bold'); axes[2].legend(fontsize=8)
    plt.tight_layout(); st.pyplot(fig); plt.close()

# ═══════════════════════════════════════════════════════════
# INTERFACE
# ═══════════════════════════════════════════════════════════
st.markdown('<div class="main-title">🧬 Analyse des Communautés de Gènes — Réseaux PPI</div>', unsafe_allow_html=True)
st.markdown('<div class="subtitle">Master 1 IES-D3S · Université Paris Nanterre · GA-PPI-Net · Soumahoro Moriba Nouho</div>', unsafe_allow_html=True)
st.divider()

# SIDEBAR
with st.sidebar:
    st.markdown("## ⚙️ Paramètres")
    st.markdown("### 📂 Chargement des fichiers")
    fichier_csv = st.file_uploader("① Réseau PPI (.csv)", type=['csv'],
                                   help="HSFinalSIVF__1_.csv — paires de gènes avec interaction et similarité")
    fichier_txt = st.file_uploader("② Communautés (.txt)", type=['txt'],
                                   help="bestSample.txt — 99 communautés détectées par GA-PPI-Net")
    st.divider()
    st.markdown("### 🎚️ Seuils")
    seuil_interaction = st.slider("Seuil interaction (∇I)", 100, 900, 150, 50,
                                   help="Score STRING minimum pour tracer un lien entre deux gènes")
    st.divider()
    st.markdown("### 📊 Mode")
    mode = st.radio("Mode d'analyse", ["🔬 Analyse d'une communauté", "⚖️ Comparaison de deux communautés"])
    st.divider()
    st.markdown("""
    **À propos**  
    Ce rapport s'appuie sur les travaux de la thèse de **Marwa Ben M'barek (2019)** — *Détection de communautés dans les grands réseaux : Application aux réseaux d'interactions de gènes* — Université Paris Nanterre.
    
    **Fitness :** `F(S) = 0.5×AVGSIM + 0.5×AVGInt/1000`
    """)

# CHARGEMENT
if not fichier_csv or not fichier_txt:
    st.info("👈 **Chargez les deux fichiers dans la barre latérale gauche pour commencer.**")
    col1, col2 = st.columns(2)
    with col1:
        st.markdown("""
        **① Réseau PPI**  
        Fichier : `HSFinalSIVF__1_.csv`  
        Contient les paires de gènes avec leur score d'interaction STRING et leur similarité sémantique GS2.
        """)
    with col2:
        st.markdown("""
        **② Communautés GA-PPI-Net**  
        Fichier : `bestSample.txt`  
        Contient les 99 communautés de gènes détectées par l'algorithme génétique GA-PPI-Net.
        """)
    st.stop()

with st.spinner("⏳ Chargement du réseau PPI..."):
    df_ppi, dict_ppi = charger_reseau(fichier_csv)

communautes = charger_communautes(fichier_txt)
st.success(f"✅ Réseau PPI : **{len(df_ppi):,} paires** chargées · **{len(communautes)} communautés** disponibles")
st.divider()

# ═══════════════════════════════════════════════════════════
# MODE 1 — ANALYSE
# ═══════════════════════════════════════════════════════════
if "Analyse" in mode:
    st.markdown("## 🔬 Analyse d'une communauté")
    num = st.selectbox("Choisir une communauté",
                       options=list(range(len(communautes))),
                       format_func=lambda i: f"Communauté {i+1}  ({len(communautes[i])} gènes)  — Gènes : {', '.join(communautes[i][:3])}...")
    comm = communautes[num]

    c1, c2, c3 = st.columns(3)
    c1.metric("Gènes", len(comm))
    c2.metric("Communauté n°", num + 1)
    c3.metric("Seuil interaction", seuil_interaction)

    with st.expander("👁️ Voir la liste complète des gènes"):
        st.write(", ".join(comm))

    if st.button("🚀 Lancer l'analyse", type="primary", use_container_width=True):
        with st.spinner("⏳ Construction du graphe et calcul des métriques..."):
            G   = construire_graphe(comm, dict_ppi, seuil_interaction)
            moy = calculer_moyennes(comm, dict_ppi)
            met = calculer_metriques(G)
            moy['degres_vals'] = list(met.get('distribution_degres', {}).values())

        fs = round(0.5 * moy['AVGSIM'] + 0.5 * moy['AVGInteraction'] / 1000, 4)

        # Mesures biologiques
        st.markdown('<div class="section-header">🧬 Mesures biologiques</div>', unsafe_allow_html=True)
        b1, b2, b3, b4 = st.columns(4)
        b1.metric("AVGSIM",          f"{moy['AVGSIM']:.4f}",           help="Similarité sémantique GS2 moyenne")
        b2.metric("AVGInteraction",  f"{moy['AVGInteraction']:.0f}",   help="Score STRING moyen (sur 1000)")
        b3.metric("Fitness F(S)",    f"{fs:.4f}",                       help="0.5×AVGSIM + 0.5×AVGInt/1000")
        b4.metric("Paires trouvées", f"{moy['paires_trouvees']} / {moy['paires_totales']}")

        # Métriques topologiques
        st.markdown('<div class="section-header">📐 Métriques topologiques</div>', unsafe_allow_html=True)
        t1, t2, t3 = st.columns(3)
        with t1:
            st.metric("① Degré moyen",     met['degre_moyen'])
            st.metric("   Hub principal",  f"{met['degre_max']} connexions → {met['hub_principal']}")
        with t2:
            st.metric("② Densité",         met['densite'])
            st.metric("③ Distance géod.",  met.get('distance_geodesique_moy', 'N/A'))
        with t3:
            st.metric("④ Diamètre",        met.get('diametre', 'N/A'))
            st.metric("⑤ Clustering",      met['clustering_moyen'])

        cm1, cm2, cm3 = st.columns(3)
        cm1.metric("⑥ Modularité",     met.get('modularite', 'N/A'))
        cm2.metric("   Sous-groupes",  met.get('nb_sous_groupes', 'N/A'))
        cm3.metric("Score Pizzuti",    met.get('score_pizzuti', 'N/A'))

        connexe = "✅ Oui" if met.get('est_connexe') else f"⚠️ Non ({met.get('nb_composantes','?')} composantes)"
        st.info(f"**Graphe connexe :** {connexe}")

        # Graphe
        st.markdown('<div class="section-header">🗺️ Graphe de la communauté</div>', unsafe_allow_html=True)
        cg1, cg2 = st.columns([3, 1])
        with cg1:
            dessiner_graphe(G, met, f"Communauté n°{num+1}")
        with cg2:
            st.markdown("**Légende**")
            st.markdown("🔴 Hub principal (degré max)")
            st.markdown("🟡→🟠 Degré croissant")
            st.markdown("**Arêtes** : épaisseur = score interaction STRING")

        # Distributions
        st.markdown('<div class="section-header">📊 Distributions</div>', unsafe_allow_html=True)
        dessiner_histogrammes(moy)

        # DAVID
        st.markdown('<div class="section-header">🔬 Validation biologique — Outil DAVID</div>', unsafe_allow_html=True)
        st.markdown("""
        Copiez la liste ci-dessous et collez-la sur **[https://davidbioinformatics.nih.gov](https://davidbioinformatics.nih.gov)**  
        → *Start Analysis* → coller les gènes → *OFFICIAL_GENE_SYMBOL* → *Homo sapiens* → *Submit List* → *Functional Annotation Chart*
        """)
        st.code("\n".join(comm), language=None)
        st.download_button("⬇️ Télécharger la liste (.txt)",
                           data="\n".join(comm),
                           file_name=f"genes_communaute_{num+1}.txt",
                           mime="text/plain")

# ═══════════════════════════════════════════════════════════
# MODE 2 — COMPARAISON
# ═══════════════════════════════════════════════════════════
else:
    st.markdown("## ⚖️ Comparaison de deux communautés")
    cs1, cs2 = st.columns(2)
    with cs1:
        numA = st.selectbox("Communauté A", list(range(len(communautes))),
                            format_func=lambda i: f"Communauté {i+1} ({len(communautes[i])} gènes)", key='A')
    with cs2:
        numB = st.selectbox("Communauté B", list(range(len(communautes))),
                            format_func=lambda i: f"Communauté {i+1} ({len(communautes[i])} gènes)", index=1, key='B')

    if st.button("🚀 Comparer", type="primary", use_container_width=True):
        commA, commB = communautes[numA], communautes[numB]
        with st.spinner("⏳ Analyse en cours..."):
            GA  = construire_graphe(commA, dict_ppi, seuil_interaction)
            moyA = calculer_moyennes(commA, dict_ppi)
            mA   = calculer_metriques(GA)
            GB   = construire_graphe(commB, dict_ppi, seuil_interaction)
            moyB = calculer_moyennes(commB, dict_ppi)
            mB   = calculer_metriques(GB)

        fsA = round(0.5*moyA['AVGSIM']+0.5*moyA['AVGInteraction']/1000, 4)
        fsB = round(0.5*moyB['AVGSIM']+0.5*moyB['AVGInteraction']/1000, 4)

        # Tableau comparatif
        st.markdown('<div class="section-header">📋 Tableau comparatif</div>', unsafe_allow_html=True)
        df_comp = pd.DataFrame({
            'Métrique': ['Nb gènes','Nb interactions','AVGSIM','AVGInteraction','Fitness F(S)',
                         '① Degré moyen','Hub principal','② Densité','③ Dist. géodésique',
                         '④ Diamètre','⑤ Clustering','⑥ Modularité','Score Pizzuti','Connexe ?'],
            f'Communauté {numA+1}': [
                mA['nb_noeuds'], mA['nb_aretes'], moyA['AVGSIM'], moyA['AVGInteraction'], fsA,
                mA['degre_moyen'], mA['hub_principal'], mA['densite'],
                mA.get('distance_geodesique_moy','N/A'), mA.get('diametre','N/A'),
                mA['clustering_moyen'], mA.get('modularite','N/A'), mA.get('score_pizzuti','N/A'),
                '✅' if mA.get('est_connexe') else '⚠️'],
            f'Communauté {numB+1}': [
                mB['nb_noeuds'], mB['nb_aretes'], moyB['AVGSIM'], moyB['AVGInteraction'], fsB,
                mB['degre_moyen'], mB['hub_principal'], mB['densite'],
                mB.get('distance_geodesique_moy','N/A'), mB.get('diametre','N/A'),
                mB['clustering_moyen'], mB.get('modularite','N/A'), mB.get('score_pizzuti','N/A'),
                '✅' if mB.get('est_connexe') else '⚠️'],
        })
        st.dataframe(df_comp, use_container_width=True, hide_index=True)

        gagnante = numA+1 if fsA >= fsB else numB+1
        st.success(f"🏆 **Meilleure fitness F(S) : Communauté {gagnante}** ({'%.4f' % max(fsA, fsB)})")

        # Graphes côte à côte
        st.markdown('<div class="section-header">🗺️ Graphes</div>', unsafe_allow_html=True)
        cg1, cg2 = st.columns(2)
        with cg1:
            st.markdown(f"**Communauté {numA+1}**")
            dessiner_graphe(GA, mA, f"Communauté {numA+1}")
        with cg2:
            st.markdown(f"**Communauté {numB+1}**")
            dessiner_graphe(GB, mB, f"Communauté {numB+1}")

        # Gènes communs
        communs = set(commA) & set(commB)
        st.markdown('<div class="section-header">🔗 Gènes en commun</div>', unsafe_allow_html=True)
        if communs:
            st.metric("Gènes partagés", len(communs))
            st.write(", ".join(sorted(communs)))
        else:
            st.info("Ces deux communautés n'ont aucun gène en commun.")
