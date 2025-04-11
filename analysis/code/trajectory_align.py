"""

"""

from PIL import Image
import glob

# 
traj_files = glob.glob("output/trajectories/*.png")

files_F = []
files_M = []
files_FM = []

for file in traj_files:
    parts = file.split("_")
    if len(parts) > 2:  # Ensure there are enough parts
        if parts[2] == "F":
            files_F.append(file)
        elif parts[2] == "M":
            files_M.append(file)
        else:
            files_FM.append(file)


# Load and resize images

def combine_image(image_list, treatment, ROWS, COLS):
  images = [Image.open(f) for f in image_list]

  # Create a blank canvas
  final_img = Image.new("RGB", (images[0].size[0]*COLS, images[0].size[1]*ROWS), "white")
  
  # Paste images into grid
  for idx, img in enumerate(images):
      row, col = divmod(idx, COLS)
      if row >= ROWS:  # Prevent overflow if extra images
          break
      x_offset = col * images[0].size[0]
      y_offset = row * images[0].size[1]
      final_img.paste(img, (x_offset, y_offset))
  
  final_img.save("output/trajectories_" + treatment + ".png", dpi=(300, 300))
  return 0


combine_image(files_F, "female", 8, 3)
combine_image(files_M, "male", 7, 3)
combine_image(files_FM,"pair", 7, 3)
