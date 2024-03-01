import tkinter as tk
from tkinter import PhotoImage, ttk, Toplevel
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


def open_tutorial():
    if not hasattr(open_tutorial, "is_open") or not open_tutorial.is_open:
        tutorial = Toplevel(app)
        tutorial.title("New Window 1")
        # This line re-enables the button when the window is closed.
        tutorial.protocol("WM_DELETE_WINDOW", lambda: on_close_tutorial(tutorial))
        open_tutorial.is_open = True
        tutorial_button["state"] = "disabled"

def open_legend():
    if not hasattr(open_legend, "is_open") or not open_legend.is_open:
        legend = Toplevel(app)
        legend.title("New Window 2")
        # This line re-enables the button when the window is closed.
        legend.protocol("WM_DELETE_WINDOW", lambda: on_close_legend(legend))
        open_legend.is_open = True
        legend_button["state"] = "disabled"


def on_close_tutorial(window):
    open_tutorial.is_open = False
    tutorial_button["state"] = "normal"
    window.destroy()

def on_close_legend(window):
    open_legend.is_open = False
    legend_button["state"] = "normal"
    window.destroy()


# Create the main application window
app = tk.Tk()
app.title("Peptide Planner")
app.geometry("800x650")
app.configure(background="white")
style = ttk.Style() # Create a style object for later use


# Pack the top_label first
top_label = tk.Label(app, background="white", text="Welcome to the Peptide Planner!", font=("Roboto", 16))
top_label.pack(pady=10)

## New window buttons
# Create a frame to hold the buttons
button_frame = tk.Frame(app)
button_frame.pack(pady=(0, 10))  # Add some padding below the frame for spacing
# Tutorial button created and packed
tutorial_button = ttk.Button(button_frame, text="Open Tutorial", command=open_tutorial)
tutorial_button.pack(side="left", fill="x", expand=True)
# Legend button created and packed
legend_button = ttk.Button(button_frame, text="Open Legend", command=open_legend)
legend_button.pack(side="right", fill="x", expand=True)

# Input text box
entry = tk.Entry(app, width=50, background="light gray", font=('Arial 12'))
entry.pack(padx=10, pady=10)

## Generate Peptide Button
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
app.bind('<Return>', on_button_click)   # enables the enter key to trigger the button


## Tabs generation
# Configure the design of the tabs
style.configure("Custom.TNotebook.Tab", font=('Helvetica', 10, 'bold'), padding=[20, 10], background='lightgray', foreground='blue')
style.map("Custom.TNotebook.Tab",
          background=[("selected", "darkblue"), ("active", "blue")], 
          foreground=[("active", "blue")])
# Create a Notebook to hold all the tabs
notebook = ttk.Notebook(app, style="Custom.TNotebook")
notebook.pack(fill='both', expand=True, padx=10, pady=10)


### Tab 1 for image and mass/charge info
tab1 = ttk.Frame(notebook)
notebook.add(tab1, text='Structure, Mass, Charge')

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


### Tab 2 for mass spectrometry
tab2 = ttk.Frame(notebook)
notebook.add(tab2, text='Mass Spectrometry')


### Tab 3 for _________________
tab3 = ttk.Frame(notebook)
notebook.add(tab3, text='Third Tab')


### Tab 4 for _________________
tab4 = ttk.Frame(notebook)
notebook.add(tab4, text='Fourth Tab')


### Tab 5 for _________________
tab5 = ttk.Frame(notebook)
notebook.add(tab5, text='Fifth Tab')


# Start the application
app.mainloop()
