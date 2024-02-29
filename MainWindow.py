import tkinter as tk
from tkinter import PhotoImage, ttk
from ImageGenerator import generate_peptide_image
from Calculations import calculate_mass, calculate_charge

# Function to be called when the button is clicked
def on_button_click(event=None):
    amino_acid_string = entry.get()
    # Generates and displays the image in the first tab
    image_path = generate_peptide_image(amino_acid_string)
    photo = PhotoImage(file=image_path)
    image_label.config(image=photo)
    image_label.image = photo  # Keep a reference
    
    # Calculates and displays the mass and net charge in the first tab
    mass = calculate_mass(amino_acid_string)
    mass_label.config(text=f"Mass: {mass}")  # Update mass label
    charge = calculate_charge(amino_acid_string)
    charge_label.config(text=f"Net Charge: {charge}")  # Update charge label


# Create the main application window
app = tk.Tk()
app.title("Peptide Planner")
app.geometry("800x600")
app.configure(background="white")


## Top text and input box
# First label on the window
top_label = tk.Label(app, background="white", text="Welcome to the Peptide Planner!", font=("Verdana", 16))
top_label.pack(pady=10)
# Input text box
entry = tk.Entry(app, width=50, background="light gray", font=('Arial 12'))
entry.pack(padx=10, pady=10)


## Button
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
# Create the button and apply the style to the button
button = ttk.Button(app, text="Generate Peptide", style="Modern.TButton", command=on_button_click)
button.pack(pady=10)
# enables the enter key to trigger the button
app.bind('<Return>', on_button_click)


## Tabs
notebook = ttk.Notebook(app)
notebook.pack(fill='both', expand=True, padx=10, pady=10)
# Tab 1 for image and mass/charge info
tab1 = ttk.Frame(notebook)
notebook.add(tab1, text='Peptide Info')
# Remaining tabs
for i in range(2, 6):       # LATER: ORGANIZE THIS FILE BY TAB
    tab = ttk.Frame(notebook)
    notebook.add(tab, text=f'Tab {i}')


## Image generation
# Frame to hold the image
image_frame = tk.Frame(tab1)
image_frame.configure(background="white")
image_frame.pack(fill='both', expand=True, padx=10, pady=10)
# Initialize the image_label with an empty image or placeholder
placeholder_image = tk.PhotoImage()  # A blank PhotoImage object
image_label = tk.Label(image_frame, image=placeholder_image)
image_label.photo = placeholder_image  # Keep a reference
image_label.pack(padx=10, pady=10)


## Mass and Net Charge
# Create a dedicated frame for the mass/charge labels
info_frame = tk.Frame(tab1)
info_frame.pack(fill='x', padx=20, pady=10)
# Mass information
mass_label = tk.Label(info_frame, text="Mass:", font=("Arial", 14), anchor="w")
mass_label.pack(side='top', pady=5, fill="x")
# Net Charge information
charge_label = tk.Label(info_frame, text="Net Charge:", font=("Arial", 14), anchor="w")
charge_label.pack(side='top', pady=5, fill="x")

# Start the application
app.mainloop()
