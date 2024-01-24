import tkinter as tk
from tkinter import ttk
from tkinter import PhotoImage
from AminoAcids import amino_acid_dict
from PyPeptBuilder import generate_peptide_image

# Function to be called when the button is clicked
def on_button_click():
    amino_acid_string = entry.get() 
    image_path = generate_peptide_image(amino_acid_string)
    photo = PhotoImage(file=image_path)
    image_label.config(image=photo)
    image_label.image = photo  # Keep a reference


# Create the main application window
app = tk.Tk()
app.title("Peptide Planner")
app.minsize(800, 600)

## Create and pack the initial label, input box, and Button
# First label on the window
top_label = tk.Label(app, text="Welcome to the Peptide Planner!")
top_label.pack(pady=10)
# Input text box
entry = tk.Entry(app, width=50, font=('Arial 12'))
entry.pack(padx=10, pady=10)  # To Do Later: restrict what user's can input
# Button
button = tk.Button(app, text="Generate Peptide", command=on_button_click)
button.pack(pady=10)


## Image generation
# Load the image (use an appropriate file path)
image_path = ""  # Replace with your image file path
photo = PhotoImage(file=image_path)

# Frame to hold the image and labels
image_frame = tk.Frame(app)
image_frame.pack(padx=10, pady=10)

# Create a label to display the image and pack it to the left
image_label = tk.Label(image_frame, image=photo)
image_label.pack(side=tk.LEFT, padx=10)

# Frame to hold the text labels
label_frame = tk.Frame(image_frame)
label_frame.pack(side=tk.LEFT, padx=20)

# Create and pack the Mass and Net Charge Labels
label1 = tk.Label(label_frame, text="Mass:   ", font=("Arial", 14), anchor="w")
label1.pack(pady=5, fill="x")
label2 = tk.Label(label_frame, text="Net Charge:   ", font=("Arial", 14), anchor="w")
label2.pack(pady=5, fill="x")


# Start the application
app.mainloop()




################################## code that might be useful later

### Combo Boxes to collect input
# # Frame to hold the first row of combo boxes
# combo_frame_one = tk.Frame(app)
# combo_frame_one.pack(pady=10)  # Adding some padding for better spacing

# # Create and pack the top 5 combo boxes
# for i in range(5):
#     # combo = ttk.Combobox(combo_frame_one, values=[f"Option {j+1}" for j in range(5)])  # past way to populate combobox
#     combo = ttk.Combobox(combo_frame_one, values=list(amino_acid_dict.keys()))
#     combo.pack(side=tk.LEFT, padx=5)  # Pack each combo box side by side with some padding


### Various Button Click stuff 
# def on_button_click():
#     # button_label.config(text="Button Clicked!")
#     # selected_value = combo.get()      # use later to display whole AA sequence
#     # label.config(text=f"Selected: {selected_value}")
#     amino_acid_string = entry.get()  # Assuming 'entry' is your input widget
