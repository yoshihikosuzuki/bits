import matplotlib.pyplot as plt
import matplotlib.image as img


def show_image(img_fname: str, width: int = 11, height: int = 11):
    """Show an image file in Jupyter Notebook."""
    fig, ax = plt.subplots(figsize=(width, height))
    ax.tick_params(labelbottom=False, bottom=False)
    ax.tick_params(labelleft=False, left=False)
    plt.imshow(img.imread(img_fname))
    plt.show()
