# 🧬 Analyse des Communautés de Gènes — Réseaux PPI & GA-PPI-Net

> **Analyse topologique et validation biologique de communautés de gènes détectées par l'algorithme génétique GA-PPI-Net dans le réseau d'interactions protéine-protéine Homo sapiens.**

[![Streamlit App](https://static.streamlit.io/badges/streamlit_badge_black_white.svg)](https://analyse-ppi-gappi-wm6zbewb7b4p4cpaxvovym.streamlit.app)
[![Python](https://img.shields.io/badge/Python-3.10+-blue.svg)](https://python.org)
[![NetworkX](https://img.shields.io/badge/NetworkX-3.x-orange.svg)](https://networkx.org)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

---

## 🔗 Application interactive en ligne

**👉 [Accéder à l'application Streamlit](https://analyse-ppi-gappi-wm6zbewb7b4p4cpaxvovym.streamlit.app)**

Aucune installation requise — chargement automatique des données au démarrage.

---

## 📌 Contexte du projet

Ce projet est réalisé dans le cadre du **stage de Master 1 IES-D3S** à l'Université Paris Nanterre (2025–2026), sous la direction de **Mme Sana Ben Hamida Mrabet**.

Il s'inscrit dans la continuité de la thèse de doctorat de **Marwa Ben M'barek (2019)** :
> *Détection de communautés dans les grands réseaux : Application aux réseaux d'interactions de gènes — Université Paris Nanterre.*

L'algorithme **GA-PPI-Net** a produit 99 communautés de gènes sur le réseau PPI humain complet (892 022 paires). Ce projet analyse, visualise et valide biologiquement ces communautés.

---

## 🎯 Objectifs

| # | Objectif |
|---|---|
| 1 | Analyser les graphes des communautés de gènes via **NetworkX** |
| 2 | Calculer 6 métriques topologiques par communauté |
| 3 | Calculer les mesures biologiques : AVGSIM, AVGInteraction, Fitness F(S) |
| 4 | Valider biologiquement 3 communautés avec **DAVID Bioinformatics** |
| 5 | Développer une interface web interactive avec **Streamlit** |

---

## 🔬 Résultats clés

| Communauté | Gènes | Hub | Fitness F(S) | Modularité | Reactome |
|---|---|---|---|---|---|
| **Comm. n°1** | 72 | GAPDH (degré 49) | 0.421 | 0.350 | 75.00% |
| **Comm. n°57** | 68 | GAPDH (degré 45) | 0.442 | 0.321 | 80.00% |
| **Comm. n°78** | 34 | MAPK14 (degré 12) | **0.446 ✅** | **0.442 ✅** | **85.29% ✅** |

**Voies biologiques de référence retrouvées :**
- ✅ **Oocyte meiosis** — Communauté n°57 (p = 0.025)
- ✅ **Apoptosis** — Communauté n°78 via Gene Ontology (p = 0.031)

---

## 🚀 Fonctionnalités de l'interface

### Mode 1 — Analyse d'une communauté
- Sélection parmi les **99 communautés** de GA-PPI-Net
- **Graphe NetworkX** coloré par degré (hub en rouge)
- Calcul automatique des **6 métriques topologiques**
- Mesures biologiques : AVGSIM, AVGInteraction, Fitness F(S)
- Histogrammes des distributions
- Export de la liste de gènes pour **DAVID**

### Mode 2 — Comparaison de deux communautés
- Tableau comparatif complet
- Deux graphes côte à côte
- Identification des gènes communs

---

## 🗂️ Structure du projet

```
analyse-ppi-gappi/
├── app.py                              # Application Streamlit
├── requirements.txt                    # Dépendances Python
├── bestSample.txt                      # 99 communautés GA-PPI-Net
└── HSFinalSIVF_communities_only.csv    # Réseau PPI filtré (95 533 interactions)
```

---

## ⚙️ Technologies utilisées

| Bibliothèque | Rôle |
|---|---|
| `streamlit` | Interface web interactive |
| `networkx` | Construction et analyse des graphes PPI |
| `matplotlib` | Visualisation des graphes et histogrammes |
| `pandas` | Manipulation du réseau PPI (CSV) |
| `numpy` | Calculs statistiques |
| `python-louvain` | Calcul de la modularité (algorithme de Louvain) |
| `scipy` | Algorithmes de disposition des graphes |

---

## 📊 Données biologiques

| Fichier | Description | Taille |
|---|---|---|
| `bestSample.txt` | 99 communautés de gènes (GA-PPI-Net) | 136 Ko |
| `HSFinalSIVF_communities_only.csv` | Réseau PPI Homo sapiens filtré | 3.1 Mo |

> Le réseau PPI complet contient **892 022 paires** de gènes humains avec pour chaque paire : le score d'interaction STRING (ValueInteraction) et le score de similarité sémantique GS2.

---

## 🧮 Métriques calculées

### Biologiques
- **AVGSIM** : similarité sémantique GS2 moyenne entre toutes les paires de gènes
- **AVGInteraction** : score STRING moyen entre toutes les paires
- **Fitness F(S)** = 0.5 × AVGSIM + 0.5 × AVGInteraction/1000

### Topologiques (NetworkX)
- **① Degré** : nombre de connexions par gène → identification des hubs
- **② Densité** : proportion de liens présents / liens possibles
- **③ Distance géodésique** : longueur moyenne du plus court chemin
- **④ Diamètre** : distance maximale dans le graphe
- **⑤ Clustering** : cohésion locale (proportion des voisins interconnectés)
- **⑥ Modularité** : structure communautaire (Louvain) — seuil significatif > 0.3

---

## 🏃 Lancer localement

```bash
# Cloner le dépôt
git clone https://github.com/moribanouho/analyse-ppi-gappi.git
cd analyse-ppi-gappi

# Installer les dépendances
pip install -r requirements.txt

# Lancer l'application
streamlit run app.py
```

---

## 📚 Références

- **Ben M'barek, M. (2019).** *Détection de communautés dans les grands réseaux : Application aux réseaux d'interactions de gènes.* Thèse de doctorat, Université Paris Nanterre.
- **Pizzuti, C. (2008).** *GA-Net: A Genetic Algorithm for Community Detection in Social Networks.*
- **STRING-db** : [https://string-db.org](https://string-db.org)
- **DAVID Bioinformatics** : [https://davidbioinformatics.nih.gov](https://davidbioinformatics.nih.gov)
- **Gene Ontology** : [https://geneontology.org](https://geneontology.org)

---

## 👤 Auteur

**Soumahoro Moriba Nouho**
Master 1 IES-D3S — Université Paris Nanterre
Stage 2025–2026 — Encadrante : Mme Sana Ben Hamida Mrabet

---

*Ce projet illustre l'application de l'algorithmique évolutionnaire et de l'analyse de graphes à la bioinformatique, avec une interface interactive accessible à des utilisateurs non programmeurs.*
