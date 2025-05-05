import requests
import pandas as pd
import networkx as nx
import plotly.graph_objects as go
import numpy as np
import igraph as ig
from plotly.offline import plot
import matplotlib.pyplot as plt
import colorsys
import os
from plotly.subplots import make_subplots
import io
from PIL import Image

# Define STRING API parameters
string_api_url = "https://string-db.org/api/json/network"
organism_id = 9606  # Human
genes = ["ADM", "CCL2", "AXL", "TLR7", "AIM2", "BATF2", "GBP1", "IL15RA", "IRF1", "TRIM21", "TTYH2"]

# Prepare API parameters
params = {
    "identifiers": "%0d".join(genes),  # Format gene list for API
    "species": organism_id,
    "required_score": 400,  # Default confidence score
    "network_type": "functional"
}

# Make API call to STRING
response = requests.get(string_api_url, params=params)
data = response.json()

# Convert JSON to dataframe
ppi_network = pd.DataFrame(data)

# Display first few interactions
print("First few interactions:")
print(ppi_network[["preferredName_A", "preferredName_B", "score"]].head())

# Create igraph for analysis
g_igraph = ig.Graph()

# Get unique genes with interactions
genes_with_interactions = list(set(ppi_network["preferredName_A"].tolist() + ppi_network["preferredName_B"].tolist()))

# Find genes without interactions
genes_without_interactions = [gene for gene in genes if gene not in genes_with_interactions]
print(f"\nGenes without interactions: {', '.join(genes_without_interactions)}")

# Add all genes as vertices (including those without interactions)
g_igraph.add_vertices(genes)

# Add edges with their scores
edges_igraph = []
edge_scores = []
for i in range(len(ppi_network)):
    nodeA = ppi_network["preferredName_A"][i]
    nodeB = ppi_network["preferredName_B"][i]
    # Get vertex indices instead of names
    idxA = g_igraph.vs.find(name=nodeA).index
    idxB = g_igraph.vs.find(name=nodeB).index
    edges_igraph.append((idxA, idxB))
    edge_scores.append(ppi_network["score"][i])

g_igraph.add_edges(edges_igraph)
g_igraph.es["weight"] = edge_scores

# Calculate node statistics
g_igraph.vs["degree"] = g_igraph.degree()

# Calculate clusters using walktrap algorithm
clusters = g_igraph.community_walktrap(weights="weight")
cluster_membership = clusters.as_clustering()
g_igraph.vs["cluster"] = [0] * g_igraph.vcount()
for i, cluster in enumerate(cluster_membership):
    for member in cluster:
        g_igraph.vs[member]["cluster"] = i

# Convert igraph to NetworkX for Plotly visualization
G = nx.Graph()

# Add nodes (including those without interactions)
for v in g_igraph.vs:
    G.add_node(v["name"], 
               degree=v["degree"],
               cluster=v["cluster"])

# Add edges
for e in g_igraph.es:
    source = g_igraph.vs[e.source]["name"]
    target = g_igraph.vs[e.target]["name"]
    G.add_edge(source, target, weight=e["weight"]/1000)

# Set node positions using Fruchterman-Reingold layout
pos = nx.spring_layout(G, seed=42, k=0.6)  # Increased k for better spacing

# Theme colors - using a light theme for white background
theme_colors = {
    'background': '#ffffff',  # White background
    'plot_bg': '#ffffff',     # White plot background
    'text': '#333333',        # Dark gray text for good contrast
    'grid': '#eeeeee',        # Light gray grid
    'node_border': '#333333', # Dark border for nodes
    'colorscale': 'Viridis',  # Keeping Viridis colorscale
    'sidebar_bg': '#f8f9fa'   # Light gray sidebar background
}

# Edge colors for different interaction strengths
edge_colors = ["#F88379", "#FF3131", "#D70040"]  # Light Red, Medium Red, Dark Red

# Create the figure with subplots (main plot only)
fig = go.Figure()

# EDGE TRACES - Create separate traces for each interaction strength
edge_low = []
edge_medium = []
edge_high = []
if G.edges():  # Only calculate if there are edges
    # Get min and max weights for normalization
    min_weight = min(G.edges[edge]['weight'] for edge in G.edges())
    max_weight = max(G.edges[edge]['weight'] for edge in G.edges())
    
    for edge in G.edges():
        x0, y0 = pos[edge[0]]
        x1, y1 = pos[edge[1]]
        
        # Normalize weight between 0 and 1 for width calculation
        weight = G.edges[edge]['weight']
        normalized_weight = (weight - min_weight) / (max_weight - min_weight) if max_weight > min_weight else 0.5
        
        # Add to appropriate list based on strength
        if normalized_weight < 0.33:
            edge_low.append((x0, y0, x1, y1, weight))
        elif normalized_weight < 0.66:
            edge_medium.append((x0, y0, x1, y1, weight))
        else:
            edge_high.append((x0, y0, x1, y1, weight))

# Create one trace for each interaction strength level
if edge_low:
    x_low, y_low = [], []
    hovertext_low = []
    for x0, y0, x1, y1, weight in edge_low:
        x_low.extend([x0, x1, None])
        y_low.extend([y0, y1, None])
        hovertext_low.append(f"Low Interaction\nScore: {weight:.3f}")
    
    low_edge_trace = go.Scatter(
        x=x_low, y=y_low,
        line=dict(width=2, color=edge_colors[0]),
        hoverinfo='text',
        hovertext=hovertext_low,
        mode='lines',
        name='Low Interaction',
        showlegend=True
    )
    fig.add_trace(low_edge_trace)

if edge_medium:
    x_medium, y_medium = [], []
    hovertext_medium = []
    for x0, y0, x1, y1, weight in edge_medium:
        x_medium.extend([x0, x1, None])
        y_medium.extend([y0, y1, None])
        hovertext_medium.append(f"Medium Interaction\nScore: {weight:.3f}")
    
    medium_edge_trace = go.Scatter(
        x=x_medium, y=y_medium,
        line=dict(width=4, color=edge_colors[1]),
        hoverinfo='text',
        hovertext=hovertext_medium,
        mode='lines',
        name='Medium Interaction',
        showlegend=True
    )
    fig.add_trace(medium_edge_trace)

if edge_high:
    x_high, y_high = [], []
    hovertext_high = []
    for x0, y0, x1, y1, weight in edge_high:
        x_high.extend([x0, x1, None])
        y_high.extend([y0, y1, None])
        hovertext_high.append(f"High Interaction\nScore: {weight:.3f}")
    
    high_edge_trace = go.Scatter(
        x=x_high, y=y_high,
        line=dict(width=6, color=edge_colors[2]),
        hoverinfo='text',
        hovertext=hovertext_high,
        mode='lines',
        name='High Interaction',
        showlegend=True
    )
    fig.add_trace(high_edge_trace)

# NODE TRACE
node_x = []
node_y = []
for node in G.nodes():
    x, y = pos[node]
    node_x.append(x)
    node_y.append(y)

# Calculate max degree for node size scaling
max_degree = max(G.nodes[node]['degree'] for node in G.nodes())
if max_degree == 0:  # Avoid division by zero
    max_degree = 1

# Prepare node data
node_sizes = []
node_colors = []
node_hover_texts = []
node_symbols = []

# Create a dictionary to map cluster numbers to vibrant colors
for node in G.nodes():
    degree = G.nodes[node]['degree']
    cluster = G.nodes[node]['cluster']
    
    # Scale size based on degree
    node_sizes.append(20 + (degree / max_degree * 30))
    
    # Use cluster as color index
    node_colors.append(cluster)
    
    # Use different symbol for non-interacting genes
    if node in genes_without_interactions:
        node_symbols.append('diamond')
    else:
        node_symbols.append('circle')
    
    # Create hover text with info
    hover_text = f"Gene: {node}\nConnections: {degree}\nCluster: {cluster+1}"
    if node in genes_without_interactions:
        hover_text += "\n(No interactions)"
    node_hover_texts.append(hover_text)

# Create node trace - INCREASED TEXT SIZE from 14 to 20
node_trace = go.Scatter(
    x=node_x, y=node_y,
    mode='markers+text',
    text=[node for node in G.nodes()],
    textposition="top center",
    textfont=dict(size=48, color=theme_colors['text'], family="Arial Black"),  # Increased size and bold font
    hoverinfo='text',
    hovertext=node_hover_texts,
    marker=dict(
        showscale=True,
        colorscale='Viridis',
        reversescale=True,
        color=node_colors,
        size=node_sizes,
        symbol=node_symbols,
        line=dict(width=2, color=theme_colors['node_border']),
        colorbar=dict(
            thickness=15,
            title='Cluster',
            xanchor='left',
            titleside='right',
            titlefont=dict(color=theme_colors['text']),
            tickfont=dict(color=theme_colors['text']),
            y=0.5,
            len=0.5
        )
    ),
    showlegend=False  # Hide from legend
)
fig.add_trace(node_trace)

# Update layout
fig.update_layout(
    titlefont=dict(size=48, color=theme_colors['text']),
    showlegend=True,
    legend=dict(
        font=dict(color=theme_colors['text']),
        x=0.95,  # Position in bottom right
        y=0.05,
        bgcolor='rgba(255,255,255,0.9)',
        bordercolor='rgba(0,0,0,0.2)',
        borderwidth=1,
        itemsizing='constant'  # Ensure uniform legend item sizes
    ),
    margin=dict(b=20, l=5, r=5, t=60),
    paper_bgcolor=theme_colors['background'],
    plot_bgcolor=theme_colors['plot_bg'],
    height=800,  # Increased height
)

# Update axes for main plot
fig.update_xaxes(showgrid=False, zeroline=False, showticklabels=False)
fig.update_yaxes(showgrid=False, zeroline=False, showticklabels=False)

# Configure the figure size for high resolution
fig.update_layout(
    width=2400,  # Increased width for higher resolution
    height=2000,  # Increased height for higher resolution
    font=dict(size=48)  # Increase font size for better readability in high-res
)

# Output the figure
fig.write_html('interactive_ppi_network.html')

# Define a function to save TIFF files using PIL without requiring Kaleido
def save_as_tiff(fig, filename, width=2400, height=2000, dpi=300):
    """
    Save a plotly figure as TIFF format using matplotlib as intermediary
    
    Parameters:
    -----------
    fig : plotly.graph_objects.Figure
        The plotly figure to save
    filename : str
        The filename to save the TIFF as
    width : int
        The width of the image in pixels
    height : int
        The height of the image in pixels
    dpi : int
        The dots per inch (resolution) of the image
    """
    # First save as a temporary HTML file
    temp_html = 'temp_plot.html'
    fig.write_html(temp_html, include_plotlyjs='cdn')
    
    # Use matplotlib to recreate a simple version of the network
    plt.figure(figsize=(width/dpi, height/dpi), dpi=dpi)
    
    # Draw edges
    for u, v in G.edges():
        x0, y0 = pos[u]
        x1, y1 = pos[v]
        weight = G.edges[u, v]['weight']
        # Determine edge width and color based on weight
        if weight < 0.45:  # Low
            edge_width = 1
            edge_color = edge_colors[0]
        elif weight < 0.6:  # Medium
            edge_width = 2
            edge_color = edge_colors[1]
        else:  # High
            edge_width = 3
            edge_color = edge_colors[2]
        plt.plot([x0, x1], [y0, y1], color=edge_color, linewidth=edge_width, alpha=0.7)
    
    # Draw nodes
    for node in G.nodes():
        x, y = pos[node]
        degree = G.nodes[node]['degree']
        cluster = G.nodes[node]['cluster']
        
        # Scale node size based on degree
        size = 100 + (degree / max_degree * 500)
        
        # Color based on cluster using viridis colormap
        color_tuple = plt.cm.viridis(cluster / max(1, len(cluster_membership)))
        
        # Use different marker for isolated nodes
        if node in genes_without_interactions:
            marker = 'D'  # Diamond
        else:
            marker = 'o'  # Circle
        
        plt.scatter(x, y, s=size, c=[color_tuple], marker=marker, 
                    edgecolors='black', linewidths=1, alpha=0.8)
        
        # Add gene name labels with increased font size
        plt.text(x, y + 0.05, node, fontsize=48, ha='center', va='center', 
                 fontweight='bold', bbox=dict(facecolor='white', alpha=0.7, boxstyle='round,pad=0.2'))
    
    # Set plot limits with some padding
    plt.xlim(min(node_x) - 0.2, max(node_x) + 0.2)
    plt.ylim(min(node_y) - 0.2, max(node_y) + 0.2)
    
    # Turn off axis
    plt.axis('off')
    
    # Save as TIFF
    plt.savefig(filename, format='tiff', dpi=dpi, bbox_inches='tight')
    plt.close()
    
    print(f"Saved TIFF image to {filename}")
    
    # Clean up temporary file
    if os.path.exists(temp_html):
        os.remove(temp_html)

# Attempt to use Kaleido for TIFF export if available
try:
    from plotly.io import write_image
    
    # Try to save directly as TIFF with high resolution using Kaleido
    try:
        # Save as TIF
        write_image(fig, 'ppi_network_highres.tif', scale=4, width=2400, height=2000)
        print("\nHigh-resolution TIF image has been saved to 'ppi_network_highres.tif'")
        
        # Save as TIFF
        write_image(fig, 'ppi_network_highres.tiff', scale=4, width=2400, height=2000)
        print("High-resolution TIFF image has been saved to 'ppi_network_highres.tiff'")
    except Exception as e:
        print(f"\nCould not save using Kaleido: {e}")
        print("Using alternative method to save TIFF files...")
        
        # Use the matplotlib-based alternative
        save_as_tiff(fig, 'ppi_network_highres.tif', width=2400, height=2000, dpi=300)
        save_as_tiff(fig, 'ppi_network_highres.tiff', width=2400, height=2000, dpi=300)
except ImportError:
    print("\nKaleido not available. Using alternative method to save TIFF files...")
    
    # Use the matplotlib-based alternative
    save_as_tiff(fig, 'ppi_network_highres.tif', width=2400, height=2000, dpi=300)
    save_as_tiff(fig, 'ppi_network_highres.tiff', width=2400, height=2000, dpi=300)

# Print some network statistics
print("\nNetwork Statistics:")
print(f"Number of nodes: {g_igraph.vcount()}")
print(f"Number of edges: {g_igraph.ecount()}")
print(f"Average degree: {sum(g_igraph.degree())/g_igraph.vcount():.2f}")
print(f"Density: {g_igraph.density():.4f}")

# Find most connected proteins
print("\nMost connected proteins:")
degree_df = pd.DataFrame({
    "Protein": g_igraph.vs["name"],
    "Degree": g_igraph.degree()
})
print(degree_df.sort_values("Degree", ascending=False).head(5))

# Print cluster information
print(f"\nNumber of clusters: {len(cluster_membership)}")
for i, cluster in enumerate(cluster_membership):
    if len(cluster) > 0:
        print(f"Cluster {i+1}: {', '.join([g_igraph.vs[member]['name'] for member in cluster])}")

print("\nInteractive visualization has been saved to 'interactive_ppi_network.html'")