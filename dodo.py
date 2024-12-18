import os
import numpy as np
import plotly.graph_objs as go


def task_build():
    delete = "rm -rf build"
    create = "mkdir build"
    make = "cd build; cmake ..; make"
    return {"actions": [delete, create, make], "verbosity": 2}


def task_send():
    create_zip = "zip -r sm24.zip ./src ./Polus"
    send = "scp sm24.zip edu-cmc-skmodel24-616-12@polus.hpc.cs.msu.ru:"
    # launch = "ssh -i ~/.ssh/id_rsa_hpc edu-cmc-skmodel24-616-12@polus.hpc.cs.msu.ru"
    delete_zip = "rm -rf sm24.zip"
    return {
        "actions": [create_zip, send, delete_zip],
    }


def task_figures():
    """Create hit maps for outputs."""

    def visualize_solutions(targets):
        """Create hit map for each file in a given directory.

        :param targets: One-element list with directory path."""
        # Step 1: Create a "figures" directory in the given path
        directory_path = targets[0]
        figures_dir = os.path.join(directory_path, "figures")
        os.makedirs(
            figures_dir, exist_ok=True
        )  # Creates the directory if it doesn't already exist

        # Step 2: Iterate over each file in the specified directory
        for filename in os.listdir(directory_path):
            file_path = os.path.join(directory_path, filename)

            # Ensure we only consider .txt files
            if (
                os.path.isfile(file_path)
                and filename.endswith(".txt")
                and filename != "CMakeCache.txt"
            ):
                # Step 3: Read the file content
                with open(file_path, "r") as file:
                    # Ignore the first line
                    lines = file.readlines()[2:]
                    # Convert remaining lines into a 2D list of floats
                    data = [list(map(float, line.split())) for line in lines]

                # Step 4: Create a NumPy array from the data
                data_array = np.array(data)

                # Step 5: Generate the 2D grid coordinates
                x = np.linspace(
                    0, 3, data_array.shape[1]
                )  # Adjust divisions if necessary
                y = np.linspace(
                    0, 3, data_array.shape[0]
                )  # Adjust divisions if necessary

                # Step 6: Create the plot
                fig = go.Figure(
                    data=go.Heatmap(
                        z=data_array,
                        x=x,  # Use the x coordinates for the columns
                        y=y,  # Use the y coordinates for the rows
                        colorscale=[[0, "darkblue"], [0.5, "yellow"], [1, "green"]],
                        colorbar=dict(title="Value"),
                        # showscale=False,
                    )
                )

                # Update layout
                fig.update_layout(
                    title="2D Heatmap for " + filename,
                    xaxis_title="X-axis",
                    yaxis_title="Y-axis",
                    xaxis=dict(nticks=10),
                    yaxis=dict(nticks=10),
                    width=800,
                    height=800,
                    # Uncomment the following line to reverse y-axis
                    # yaxis_autorange='reversed'
                )

                # Step 7: Save the plot as an image in the "figures" directory
                image_file_path = os.path.join(
                    figures_dir, f"{os.path.splitext(filename)[0]}.png"
                )
                fig.write_image(image_file_path)  # Save as PNG
                print(f"Saved heatmap to {image_file_path}")

    return {
        "actions": [visualize_solutions],
        "targets": ["build"],
    }
