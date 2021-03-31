import matplotlib.pyplot as plt
import matplotlib.image as img


def show_image(img_fname: str,
               size: int = 11):
    """Show an image file using matplotlib.

    positional arguments:
      @ img_fname : File name of an image.

    optional arguments:
      @ width  : In inch.
      @ height : In inch.
    """
    _, ax = plt.subplots(figsize=(size, size))
    ax.tick_params(labelbottom=False, bottom=False)
    ax.tick_params(labelleft=False, left=False)
    plt.imshow(img.imread(img_fname))
    plt.show()
