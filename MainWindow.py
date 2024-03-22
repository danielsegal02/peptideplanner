import tkinter as tk
from tkinter import PhotoImage, ttk, Toplevel
import pandas as pd
from ImageGenerator import generate_peptide_image, generate_mass_spec
from Calculations import calculate_mass, calculate_charge

# Function to be called when the button is clicked
def on_button_click(event=None):
    amino_acid_string = entry.get()
    # Generates and displays the image in the first tab
    chem_struct_image_path = generate_peptide_image(amino_acid_string)
    chem_struct_photo = PhotoImage(file=chem_struct_image_path)
    chem_struct_image_label.config(image=chem_struct_photo)
    chem_struct_image_label.image = chem_struct_photo  # Keep a reference
    
    # Calculates and displays the mass and net charge in the first tab
    mass = calculate_mass(amino_acid_string)
    mass_label.config(text=f"Mass: {mass}")  # Update mass label
    charge = calculate_charge(amino_acid_string)
    charge_label.config(text=f"Net Charge: {charge}")  # Update charge label

    # Calculates and displays the mass spec in the second tab
    image_path = generate_mass_spec(mass, charge)
    mass_spec_photo = PhotoImage(file=image_path)
    mass_spec_image_label.config(image=mass_spec_photo)
    mass_spec_image_label.image = mass_spec_photo  # Keep a reference


def open_tutorial():
    if not hasattr(open_tutorial, "is_open") or not open_tutorial.is_open:
        tutorial = Toplevel(app)
        tutorial.title("Tutorial")
        tutorial.minsize(300, 600)
        # This line re-enables the button when the window is closed.
        tutorial.protocol("WM_DELETE_WINDOW", lambda: on_close_tutorial(tutorial))
        open_tutorial.is_open = True
        tutorial_button["state"] = "disabled"

def open_legend():
    if not hasattr(open_legend, "is_open") or not open_legend.is_open:
        legend = Toplevel(app)
        legend.title("Amino Acid Legend")
        legend.minsize(300, 600)
        
        # Create the treeview widget
        columns = ("Code", "Name")
        tree = ttk.Treeview(legend, columns=columns, show="headings")
        tree.heading("Code", text="Input Code")
        tree.heading("Name", text="Name")
        
        # Adjust the columns' width to the content
        tree.column("Code", minwidth=50, width=100, anchor='center')
        tree.column("Name", minwidth=100, width=150, anchor='center')
        
        # Read the CSV file using pandas and populate the treeview with the amino acids columns
        df = pd.read_csv("AminoAcidTable.csv")
        for index, row in df.iterrows():
            tree.insert("", tk.END, values=(row["Code"], row["Name"]))
                
        tree.pack(expand=True, fill='both')

        # Define the Add Amino Acid button with a command that includes itself
        add_AA_button = ttk.Button(legend, text="Add a new Amino Acid", command=lambda: open_add_AA_window(add_AA_button, legend))
        add_AA_button.pack(pady=10)
        
        legend.protocol("WM_DELETE_WINDOW", lambda: on_close_legend(legend))
        open_legend.is_open = True
        open_legend.add_AA_window = None  # Initially, there's no add_AA_window
        legend_button["state"] = "disabled"

def open_add_AA_window(button, legend):
    # Check if the add_AA_window is already open using a global variable or an attribute of the legend window
    if not hasattr(legend, "add_AA_window_open") or not legend.add_AA_window_open:
        add_AA_window = Toplevel(app)
        add_AA_window.title("Add Amino Acids")
        add_AA_window.geometry("300x200")
        
        label = tk.Label(add_AA_window, text="Add your amino acid details here.")
        label.pack(pady=10)

        button["state"] = "disabled"
        legend.add_AA_window_open = True  # Mark the add_AA_window as open

        # Instead of using open_add_AA_window.is_open, store the add_AA_window reference in the legend window for access during closure
        legend.add_AA_window_ref = add_AA_window

        add_AA_window.protocol("WM_DELETE_WINDOW", lambda: on_close_add_AA(add_AA_window, button, legend))

def on_close_tutorial(window):
    open_tutorial.is_open = False
    tutorial_button["state"] = "normal"
    window.destroy()

def on_close_legend(window):
    # Check if the add_AA_window is open and close it if necessary
    if hasattr(window, "add_AA_window_open") and window.add_AA_window_open:
        if hasattr(window, "add_AA_window_ref"):
            window.add_AA_window_ref.destroy()  # Close the add_AA_window if it's open
            window.add_AA_window_ref = None
    window.add_AA_window_open = False  # Reset the state
    # Reset the legend is_open state and enable the legend button
    open_legend.is_open = False
    legend_button["state"] = "normal"
    window.destroy()

def on_close_add_AA(window, button, legend):
    # Mark the add_AA_window as closed and enable the button for reopening it
    legend.add_AA_window_open = False
    button["state"] = "normal"
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


## Entry box and Termini select boxes
# Frame that entry and combo boxes will be packed in
combo_entry_frame = tk.Frame(app)
combo_entry_frame.pack(padx=10, pady=10)
# Left Combo Box
combo_options_left = ['Option 1', 'Option 2', 'Option 3']
combo_box_left = ttk.Combobox(combo_entry_frame, values=combo_options_left, state="readonly", width=15)
combo_box_left.pack(side='left', padx=(0, 5))  # Pack to the left side with some padding
# Input text box
entry = tk.Entry(combo_entry_frame, width=50, background="light gray", font=('Arial 12'))
entry.pack(side='left', padx=5)
# Right Combo Box
combo_options_right = ['A', 'B', 'C']
combo_box_right = ttk.Combobox(combo_entry_frame, values=combo_options_right, state="readonly", width=15)
combo_box_right.pack(side='left', padx=(5, 0))  # Pack to the left side, which effectively places it to the right of the entry


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
notebook.add(tab1, text='Mass & Charge')

## Chemical structure image generation
# Frame to hold the image
chem_struct_image_frame = tk.Frame(tab1)
chem_struct_image_frame.configure(background="white")
chem_struct_image_frame.pack(fill='both', expand=True, padx=10, pady=10)
# Initialize the chem_struct_image_label with an empty image or placeholder
chem_struct_placeholder_image = tk.PhotoImage()  # A blank PhotoImage object
chem_struct_image_label = tk.Label(chem_struct_image_frame, image=chem_struct_placeholder_image)
chem_struct_image_label.photo = chem_struct_placeholder_image  # Keep a reference
chem_struct_image_label.pack(padx=10, pady=10)

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

## Mass Spec plot image generation
# Frame to hold the image
mass_spec_image_frame = tk.Frame(tab2)
mass_spec_image_frame.configure(background="white")
mass_spec_image_frame.pack(fill='both', expand=True, padx=10, pady=10)
# Initialize the mass_spec_image_label with an empty image or placeholder
mass_spec_placeholder_image = tk.PhotoImage()  # A blank PhotoImage object
mass_spec_image_label = tk.Label(mass_spec_image_frame, image=mass_spec_placeholder_image)
mass_spec_image_label.photo = mass_spec_placeholder_image  # Keep a reference
mass_spec_image_label.pack(padx=10, pady=10)


### Tab 3 for _________________
tab3 = ttk.Frame(notebook)
notebook.add(tab3, text='Conjugation')


### Tab 4 for _________________
tab4 = ttk.Frame(notebook)
notebook.add(tab4, text='Structure')


# Start the application
app.mainloop()
