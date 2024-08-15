# Function to visualize the joint position and whole body COM location in the x-axis

import matplotlib.pyplot as plt
from PIL import Image

def visualizeposture(positions_x, positions_y, bodycom_x, time, title=None, plotfile=None):

    plt.cla()
    axes = plt.gca()
    axes.set_xlim([-80, 180])
    axes.set_ylim([-180, 80])

    plt.plot([positions_x[2], positions_x[6]], [positions_y[2], positions_y[6]], color='gray')
    plt.plot([positions_x[6], positions_x[7]], [positions_y[6], positions_y[7]], color='gray')
    plt.plot([positions_x[7], positions_x[8]], [positions_y[7], positions_y[8]], color='gray')

    plt.plot([0, positions_x[0]], [0, positions_y[0]], color='black')
    plt.plot([positions_x[0], positions_x[1]], [positions_y[0], positions_y[1]], color='black')
    plt.plot([positions_x[1], positions_x[2]], [positions_y[1], positions_y[2]], color='black')
    plt.plot([positions_x[2], positions_x[3]], [positions_y[2], positions_y[3]], color='black')
    plt.plot([positions_x[3], positions_x[4]], [positions_y[3], positions_y[4]], color='black')
    plt.plot([positions_x[4], positions_x[5]], [positions_y[4], positions_y[5]], color='black')
    plt.scatter(bodycom_x, positions_y[5], color='green')
    plt.scatter(0, 0)

    plt.xlabel("X-Position (cm)", fontsize=15)
    plt.ylabel("Y-Position (cm)", fontsize=15)
    if title is not None:
        plt.title(title, fontsize=15, fontweight="bold", loc='left')
    plt.xticks(fontsize=11)
    plt.yticks(fontsize=11)

    if plotfile is not None:
        plt.savefig(plotfile)
    plt.pause(time)

    return 'Visualization'


def combine_images(columns, space, images):
    # Code taken from here: https://stackoverflow.com/questions/72723928/how-to-combine-several-images-to-one-image-in-a-grid-structure-in-python

    rows = len(images) // columns
    if len(images) % columns:
        rows += 1
    width_max = max([Image.open(image).width for image in images])
    height_max = max([Image.open(image).height for image in images])
    background_width = width_max*columns + (space*columns)-space
    background_height = height_max*rows + (space*rows)-space
    background = Image.new('RGBA', (background_width, background_height), (255, 255, 255, 255))
    x = 0
    y = 0
    for i, image in enumerate(images):
        img = Image.open(image)
        x_offset = int((width_max-img.width)/2)
        y_offset = int((height_max-img.height)/2)
        background.paste(img, (x+x_offset, y+y_offset))
        x += width_max + space
        if (i+1) % columns == 0:
            y += height_max + space
            x = 0
    background.save('Figures/Fig1-PosturePrediction/Figure1.png')
