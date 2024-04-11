import tkinter as tk
from tkinter import PhotoImage, ttk, Toplevel
import pandas as pd
from ImageGenerator import generate_peptide_image, generate_mass_spec
from Calculations import calculate_mass, calculate_charge, calculate_reagent_mass, calculate_reagent_volume, calculate_solvent_volume

# Function to be called when the button is clicked
def handle_generate_button_click(event=None):
    ## Clear any previous outputs first
    # Reset images to placeholder or clear them
    chem_struct_image_label.config(image=chem_struct_placeholder_image)
    mass_spec_image_label.config(image=mass_spec_placeholder_image)
    # Reset text labels for mass and charge to default or empty values
    mass_label.config(text="Mass:")
    charge_label.config(text="Net Charge for MS:")

    # Pulls user inputs
    amino_acid_string = main_entry_box.get()
    n_terminus = combo_box_left.get()
    c_terminus = combo_box_right.get()

    ## Error catching
    # Removes trailing spaces from input
    amino_acid_string = amino_acid_string.rstrip(" ")
    # Clear the entry box
    main_entry_box.delete(0, tk.END)
    # Insert the input without the trailing spaces into the entry box
    main_entry_box.insert(0, amino_acid_string)

    # Checks if the input is empty
    if amino_acid_string == "":
        error_msg.config(text="Error: Bad Input! Please check legend, exclude spaces, and ensure input is valid.")
        return

    # Checks if there are any characters in the input that are not in the csv
    amino_acid_df = pd.read_csv("AminoAcidTable.csv")   # Load amino acid data from a CSV file into a DataFrame
    amino_acid_codes = amino_acid_df["Code"].tolist()       # Convert the 'Code' column of the DataFrame to a list for easy lookup
    # Iterate over each character in the given amino acid string
    for amino_acid in amino_acid_string:
        # Check if the current amino acid code is not in the list of valid codes
        if amino_acid not in amino_acid_codes:
            # If an invalid code is found, display an error message and exit the function
            error_msg.config(text="Error: Bad Input! Please check legend, exclude spaces, and ensure input is valid.")
            return

    # If no errors are found, clear the error message
    error_msg.config(text="")

    # Generates and displays the image in the first tab
    chem_struct_image_path = generate_peptide_image(amino_acid_string, n_terminus, c_terminus)
    chem_struct_photo = PhotoImage(file=chem_struct_image_path)
    chem_struct_image_label.config(image=chem_struct_photo)
    chem_struct_image_label.image = chem_struct_photo  # Keep a reference
    
    # Calculates and displays the mass and net charge in the first tab
    mass = calculate_mass(amino_acid_string, n_terminus, c_terminus)
    mass_label.config(text=f"Mass: {mass}")  # Update mass label
    charge = calculate_charge(amino_acid_string, n_terminus, c_terminus)
    charge_label.config(text=f"Net Charge for MS: {charge}")  # Update charge label

    # Calculates and displays the mass spec in the second tab
    image_path = generate_mass_spec(mass, charge)
    mass_spec_photo = PhotoImage(file=image_path)
    mass_spec_image_label.config(image=mass_spec_photo)
    mass_spec_image_label.image = mass_spec_photo  # Keep a reference


def open_tutorial():
    if not hasattr(open_tutorial, "is_open") or not open_tutorial.is_open:
        tutorial = Toplevel(app)
        tutorial.title("Tutorial")
        tutorial.geometry("600x600")
        tutorial.protocol("WM_DELETE_WINDOW", lambda: on_close_tutorial(tutorial))
        open_tutorial.is_open = True
        tutorial_button["state"] = "disabled"

        # Tutorial Title
        tutorial_title = tk.Label(tutorial, text="Tutorial", font=("Arial", 18, "bold"))
        tutorial_title.pack(pady=(10, 5))

        # Tutorial Content
        tutorial_text = """
1. In the gray text box toward the top of the screen, enter the desired amino acid sequence. Be sure to type each amino acid using its corresponding letter/symbol as listed in the Legend. (For reference, you can click “Open Legend”)

2. Click “Generate Peptide” to populate the functions under all tabs. The chemical structure, mass, and net charge of the peptide should appear on the first page of the application.

3. To view the predicted mass spectrometry graph, click the “Mass Spectrometry” tab.

4. To do reagent calculations, click the “Conjugation” tab.

5. To view the predicted secondary structure, click the “Structure” tab.

6. To enter a new amino acid sequence, clear the gray text box at the top and type in the new sequence. Click “Generate Peptide”.
        """

        # Create a Text widget for the tutorial content
        tutorial_content = tk.Text(tutorial, wrap="word", font=("Arial", 12), bg="light grey", borderwidth=2, relief="solid")
        tutorial_content.insert("1.0", tutorial_text)

        # Define a tag for bold and slightly larger font for the numbers
        tutorial_content.tag_configure("list_number", font=("Arial", 14, "bold"))

        # Find all occurrences of a pattern (in this case, the numbers followed by a period)
        import re
        for match in re.finditer(r"\d+\.", tutorial_text):
            start = "1.0 + {} chars".format(match.start())
            end = "1.0 + {} chars".format(match.end())
            tutorial_content.tag_add("list_number", start, end)

        tutorial_content.config(state="disabled")  # Make the text widget read-only
        tutorial_content.pack(padx=10, pady=10, fill="both", expand=True)

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
app.geometry("800x700")
app.configure(background="white")
style = ttk.Style() # Create a style object for later use


# Pack the top_label first
top_label = tk.Label(app, background="white", text="Welcome to the Peptide Planner!", font=("Roboto", 16))
top_label.pack(pady=10)

## New window buttons
# Create a frame to hold the buttons
button_frame = tk.Frame(app)
button_frame.pack(pady=(0, 10))  # Add some padding below the frame for spacing
style.configure("TButton", font=('Helvetica', 12), width=20) # Creates a style element for the top two buttons
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
combo_options_left = ['Amine', 'Acetyl']
combo_box_left = ttk.Combobox(combo_entry_frame, values=combo_options_left, state="readonly", width=15)
combo_box_left.current(0)
combo_box_left.pack(side='left', padx=(0, 5))  # Pack to the left side with some padding
# Input text box
main_entry_box = ttk.Entry(combo_entry_frame, width=50, background="light gray", font=('Arial 12'))
main_entry_box.pack(side='left', padx=5)
# Right Combo Box
combo_options_right = ['Amide', 'Acid']
combo_box_right = ttk.Combobox(combo_entry_frame, values=combo_options_right, state="readonly", width=15)
combo_box_right.current(0)
combo_box_right.pack(side='left', padx=(5, 0))  # Pack to the left side, which effectively places it to the right of the entry


# Error Message Label
error_msg = tk.Label(app, background="white", foreground= 'red', text=" ", font=("Roboto", 10))
error_msg.pack(pady=5)


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
generate_button = ttk.Button(app, text="Generate Peptide", style="Modern.TButton", command=handle_generate_button_click)
generate_button.pack(pady=10)
app.bind('<Return>', handle_generate_button_click)   # enables the enter key to trigger the button


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
charge_label = tk.Label(info_frame, text="Net Charge for MS:", font=("Arial", 14), anchor="w")
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


### Tab 3 for Reagents
tab3 = ttk.Frame(notebook)
notebook.add(tab3, text='Conjugation')

# Add a label at the top of tab3 for tutorial reference
help_label = ttk.Label(tab3, text="*See tutorial for any help calculating the values*", foreground="blue", font=('Roboto', 11), anchor="center")
help_label.pack(padx=10, pady=(10, 0))

# Adjust style elements for a nicer look
style = ttk.Style()
style.configure("TLabel", font=('Helvetica', 12), anchor="center")
style.configure("TEntry", font=('Helvetica', 12))
style.configure("TRadiobutton", font=('Helvetica', 12))

# Frame for central column in tab3 with style adjustments
reagents_frame = ttk.Frame(tab3)
reagents_frame.pack(padx=10, pady=10, expand=True)
reagents_frame.grid_columnconfigure(0, weight=1) # Configure the frame to use all available space and center contents

## Top answer labels
# Define answer label text specifics
conjugation_answers = [
    {"text_label": "Reagent Mass", "answer_text": "_________________", "units": "mg", "key": "reagent_mass"},
    {"text_label": "Reagent Volume", "answer_text": "_________________", "units": "μL", "key": "reagent_vol"},
    {"text_label": "Solvent Volume", "answer_text": "_________________", "units": "mL", "key": "solvent_vol"},
]
# Initialize dictionaries to store label references
text_labels = {}
answer_labels = {}
units_labels = {}
# Iterate over the specifications to create, grid, and store elements
for index, spec in enumerate(conjugation_answers):
    # Text Label
    text_labels[spec["key"]] = ttk.Label(reagents_frame, text=spec["text_label"], style="TLabel")
    text_labels[spec["key"]].grid(column=0, row=index, sticky="EW")

    # Answer Label
    answer_labels[spec["key"]] = ttk.Label(reagents_frame, text=spec["answer_text"], style="TLabel")
    answer_labels[spec["key"]].grid(column=1, row=index, sticky="EW")

    # Units Label
    units_labels[spec["key"]] = ttk.Label(reagents_frame, text=spec["units"], style="TLabel")
    units_labels[spec["key"]].grid(column=2, row=index, sticky="EW")
# Now, you can access a specific label like so:
# text_labels['reagent_mass'].configure(text="New Reagent Mass")  # To update the text of the reagent mass text label

## Radio Buttons
# Creates a subframe for the radio buttons
radio_frame = ttk.Frame(reagents_frame)
radio_frame.grid(column=1, row=3, sticky="EW", pady=(15, 5))  # Adds 10px padding above the radio frame
radio_frame.grid_columnconfigure(0, weight=1)
radio_frame.grid_columnconfigure(1, weight=1)

# Creating the radio buttons for Dry or Wet
dry_wet_var = tk.StringVar(value="Dry")  # Sets the default to "Dry"
ttk.Radiobutton(radio_frame, text="Dry", variable=dry_wet_var, value="Dry", style="TRadiobutton").grid(column=0, row=0)
ttk.Radiobutton(radio_frame, text="Wet", variable=dry_wet_var, value="Wet", style="TRadiobutton").grid(column=1, row=0)

# Define the callback function for enabling/disabling the "Reagent density" entry box with the radio buttons
def toggle_reagent_density(*args):
    if dry_wet_var.get() == "Dry":
        entry_boxes['reagent_density'].configure(state="disabled", background="light grey")
    else:
        entry_boxes['reagent_density'].configure(state="normal", background="white")

# Bind the callback function to the selection of the radio button
dry_wet_var.trace_add("write", toggle_reagent_density)

## User labels and entry input boxes
# Define GUI element specifications
conjugation_entries = [
    {"name": "1. Peptide scale", "units": "mmol", "key": "pep_scale"},
    {"name": "2. % resin used", "units": "%", "key": "resin_used"},
    {"name": "3. Reagent MW", "units": "g/mol", "key": "reagent_MW"},
    {"name": "4. Reagent equiv.", "units": None, "key": "reagent_equiv"},
    {"name": "5. Reagent density", "units": "g/mL", "key": "reagent_density"},
    {"name": "6. Solvent Factor", "units": None, "key": "solvent_factor"},
]

# Initialize dictionaries to store widget references
name_labels = {}
entry_boxes = {}
units_labels = {}

# Iterate over the specifications to create, grid, and store elements
for index, elem in enumerate(conjugation_entries, start=4):  # Starting at row 4 in the center column
    # Label
    name_labels[elem["key"]] = ttk.Label(reagents_frame, text=elem["name"], style="TLabel")
    name_labels[elem["key"]].grid(column=0, row=index, sticky="W", padx=5, pady=(10, 0))

    # Entry Box
    entry_boxes[elem["key"]] = ttk.Entry(reagents_frame, font=('Helvetica', 12), width=20)
    entry_boxes[elem["key"]].grid(column=1, row=index, sticky="EW", padx=(2, 10), pady=2)

    # Units Label (if applicable)
    if elem["units"]:
        units_key = f"{elem['key']}_units"  # Create a unique key for units labels
        units_labels[units_key] = ttk.Label(reagents_frame, text=elem["units"], style="TLabel")
        units_labels[units_key].grid(column=2, row=index, sticky="W", padx=5, pady=(10, 0))

# Now, to access a specific widget, you call it like this:
# entry_boxes['pep_scale'].get()  # To get the value from the peptide scale entry box

# Manually invoke the function to apply the effect at startup
toggle_reagent_density()

# Create a function to handle the reagent calc button click
def handle_calc_button_click():
    """
    To calculate:
    - Reagent Mass:
        1. Dry must be chosen.
        2. The first four entry boxes in the tab must be filled out.
    - Reagent Volume:
        1. Wet must be chosen.
        2. The first five entry boxes in the tab must be filled out.
    - Solvent Factor: 
        1. At least the first four entry boxes in the tab must be filled out.
        Optional: It can be the first five entry boxes as well.
        2. The solvent factor entry box must be filled out.
    """
    # Reset the labels 
    answer_labels["reagent_mass"].config(text=f"_________________")
    answer_labels["reagent_vol"].config(text=f"_________________")
    answer_labels["solvent_vol"].config(text=f"_________________")


    # Extract values from entry boxes (these should always be filled out)
    pep_scale = float(entry_boxes['pep_scale'].get())
    resin_used = float(entry_boxes['resin_used'].get())
    reagent_MW = float(entry_boxes['reagent_MW'].get())
    reagent_equiv = float(entry_boxes['reagent_equiv'].get())

    # All formulas call for the persent resin used, so convert it before passing it to the functions
    percent_resin_used = resin_used / 100

    # Dry is chosen AND the first four boxes filled out --> solve for reagent mass
    if dry_wet_var.get() == "Dry":
        # Perform the calculation
        reagent_mass = calculate_reagent_mass(pep_scale, percent_resin_used, reagent_MW, reagent_equiv)
        # Update the label with the result
        answer_labels["reagent_mass"].config(text=f"{reagent_mass}")
        
    # Wet is chosen AND the first five boxes filled out--> solve for reagent volume
    if dry_wet_var.get() == "Wet":
        # Extract value the extra value needed for the calculation
        reagent_density = float(entry_boxes['reagent_density'].get())
        # Perform the calculation
        reagent_vol = calculate_reagent_volume(pep_scale, percent_resin_used, reagent_MW, reagent_equiv, reagent_density)
        # Update the label with the result
        answer_labels["reagent_vol"].config(text=f"{reagent_vol}")

    # Solvent factor filled out AND at least the first four boxes filled out--> solve for solvent volume
    if entry_boxes['solvent_factor'].get():
        # Extract value the extra value needed for the calculation
        solvent_factor = float(entry_boxes['solvent_factor'].get())
        # Perform the calculation
        solvent_vol = calculate_solvent_volume(pep_scale, percent_resin_used, solvent_factor)
        # Update the label with the result
        answer_labels["solvent_vol"].config(text=f"{solvent_vol}")

# Creating the reagent calculate button
ttk.Button(reagents_frame, text="Calculate", command=handle_calc_button_click, style="Modern.TButton").grid(column=1, row=10, sticky="EW", padx=(2, 10), pady=2)


### Tab 4 for Secondary Structure
tab4 = ttk.Frame(notebook)
notebook.add(tab4, text='Structure')

### Tab 5 for Credits
tab5 = ttk.Frame(notebook)
notebook.add(tab5, text='Credits')

# Add a Text widget to tab5 for displaying the credits information
credits_info = tk.Text(tab5, wrap="word", font=("Cambria", 12), bg="white", borderwidth=2, relief="solid", height=10)
credits_info_text = """
This application was a result of a Capstone project for Medina Labs. 

Penn State University Main Campus
2024 Spring

Group Members:
Jubrial Saigh
Arianna Parisi
Chase Burkhart
Daniel Segal
Danny Takacs
Nick Czwartacki

Github link: https://github.com/danielsegal02/peptideplanner.git\n
"""
# Inserts the credit text into the widget
credits_info.insert("end", credits_info_text)

# Inserts the group pic into the widget
image_path = "Images\dog temp image.gif"  # Replace this with the path to your image file
group_photo = PhotoImage(file=image_path)
credits_info.image_create("end", image=group_photo)  # This inserts the image at the end of the text

# hehe
easter_egg = """
\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n straight up homie
"""
credits_info.insert("end", easter_egg)

# Centers the text
credits_info.tag_configure("center", justify='center') # Create a tag to center-align text
credits_info.tag_add("center", "1.0", "end") # Apply the 'center' tag to the entire content of the widget

# Make the text widget read-only to prevent user editing
credits_info.config(state="disabled")

# Pack the Text widget in tab5
credits_info.pack(padx=10, pady=10, fill="both", expand=True)  


# Start the application
app.mainloop()
