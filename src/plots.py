import numpy as np
import plotly.graph_objects as go


def visualize_solution(file_path):
    # Step 1: Read the file content
    with open(file_path, "r") as file:
        # Ignore the first line

        lines = file.readlines()[2:]
        # Convert remaining lines into a 2D list of floats
        data = [list(map(float, line.split())) for line in lines]
    # Step 2: Create a NumPy array from the data
    data_array = np.array(data)

    # Step 3: Generate the 2D grid coordinates
    x = np.linspace(0, 3, data_array.shape[1])  # 10 divisions in the x direction
    y = np.linspace(0, 3, data_array.shape[0])  # 10 divisions in the y direction
    X, Y = np.meshgrid(x, y)

    # Step 4: Create the plot
    fig = go.Figure(
        data=go.Heatmap(
            z=data_array,
            x=x,  # Use the x coordinates for the columns
            y=y,  # Use the y coordinates for the rows
            colorscale=[[0, "darkblue"], [0.5, "yellow"], [1, "green"]],
            colorbar=dict(title="Value"),
        )
    )

    # Update layout
    fig.update_layout(
        title="2D Heatmap for " + file_path,
        xaxis_title="X-axis",
        yaxis_title="Y-axis",
        xaxis=dict(nticks=10),
        yaxis=dict(nticks=10),
        width=800,
        height=800,
        # autosize='reversed'  # Reverse y-axis to match the specified coordinates
    )
    return fig
