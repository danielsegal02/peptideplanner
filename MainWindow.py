import tkinter as tk
from tkinter import ttk
from tkinter import PhotoImage
from ImageGenerator import generate_peptide_image
from Calculations import calculate_mass, calculate_charge

# Function to be called when the button is clicked
def on_button_click(event=None):
    amino_acid_string = entry.get() 
    uppercase_AA_string = amino_acid_string.upper()
    # Update the entry with uppercase text
    entry.delete(0, tk.END)  # Remove the current text
    entry.insert(0, uppercase_AA_string)  # Insert the uppercase text

    # Generates and displays the image
    image_path = generate_peptide_image(uppercase_AA_string)
    photo = PhotoImage(file=image_path)
    image_label.config(image=photo)
    image_label.image = photo  # Keep a reference
    
    # Calculates and displays the mass and net charge
    mass = calculate_mass(uppercase_AA_string)
    mass_label.config(text=f"Mass: {mass}")  # Update mass label
    charge = calculate_charge(uppercase_AA_string)
    charge_label.config(text=f"Net Charge: {charge}")  # Update charge label


# Create the main application window
app = tk.Tk()
app.title("Peptide Planner")
app.minsize(800, 500)
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
mass_label = tk.Label(info_frame, background="white", text="Mass:", font=("Arial", 14), anchor="w")
mass_label.pack(side='top', pady=5, fill="x")
# Net Charge Info
charge_label = tk.Label(info_frame, background="white", text="Net Charge:", font=("Arial", 14), anchor="w")
charge_label.pack(side='top', pady=5, fill="x")


# Start the application
app.mainloop()
