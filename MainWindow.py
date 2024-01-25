import tkinter as tk
from tkinter import ttk
from tkinter import PhotoImage
from AminoAcids import amino_acid_dict
from PyPeptBuilder import generate_peptide_image

# Function to be called when the button is clicked
def on_button_click(event=None):
    amino_acid_string = entry.get() 
    image_path = generate_peptide_image(amino_acid_string)
    photo = PhotoImage(file=image_path)
    image_label.config(image=photo)
    image_label.image = photo  # Keep a reference
    # Here you can also update label1 and label2 with the Mass and Net Charge information


# Create the main application window
app = tk.Tk()
app.title("Peptide Planner")
app.minsize(800, 500)
app.configure(background="white")

## Create and pack the initial label, input box, and Button
# First label on the window
top_label = tk.Label(app, background="white", text="Welcome to the Peptide Planner!", font=("Verdana", 16))
top_label.pack(pady=10)
# Input text box
entry = tk.Entry(app, width=50, background="light gray", font=('Arial 12'))
entry.pack(padx=10, pady=10)  # To Do Later: restrict what user can input (i.e. empty string, format it then pass it in, etc...)
# Button
# Create a style object
style = ttk.Style()
# Configure a more modern-looking button style
style.configure("Modern.TButton",
                font=('Calibri', 12),
                background='black',
                foreground='blue',
                borderwidth=1,
                relief="raised",
                padding=3)
# Apply a layout to make the button appear with 'rounded corners'
style.layout("Modern.TButton",
             [('Button.border', {'sticky': 'nswe', 'children':
                 [('Button.padding', {'sticky': 'nswe', 'children':
                     [('Button.label', {'sticky': 'nswe'})],
                 })],
             })])
# Apply the style to the button
button = ttk.Button(app, text="Generate Peptide", style="Modern.TButton", command=on_button_click)
button.pack(pady=10)
# enables the enter key to trigger the button
app.bind('<Return>', on_button_click)


## Image generation
# Frame to hold the image
image_frame = tk.Frame(app)
image_frame.configure(background="white")
image_frame.pack(fill='both', expand=True, padx=10, pady=10)
# Initialize the image_label with an empty image or placeholder
placeholder_image = tk.PhotoImage()  # A blank PhotoImage object
image_label = tk.Label(image_frame, image=placeholder_image)
image_label.photo = placeholder_image  # Keep a reference
image_label.pack(padx=10, pady=10)


## Create a dedicated frame for the Mass and Net Charge Labels
info_frame = tk.Frame(app)
info_frame.configure(background="white")
info_frame.pack(fill='x', padx=20, pady=10)  # This will make sure it expands with the window
# Mass Information
label1 = tk.Label(info_frame, background="white", text="Mass:   ", font=("Arial", 14), anchor="w")
label1.pack(side='top', pady=5, fill="x")
# Net Charge Info
label2 = tk.Label(info_frame, background="white", text="Net Charge:   ", font=("Arial", 14), anchor="w")
label2.pack(side='top', pady=5, fill="x")


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
