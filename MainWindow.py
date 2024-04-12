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
        tutorial.geometry("700x700")
        tutorial.protocol("WM_DELETE_WINDOW", lambda: on_close_tutorial(tutorial))
        open_tutorial.is_open = True
        tutorial_button["state"] = "disabled"

        # Tutorial Title
        tutorial_title = tk.Label(tutorial, text="Tutorial", font=("Arial", 18, "bold"))
        tutorial_title.pack(pady=(10, 5))

        # Tutorial Content
        tutorial_text = """
NOTE: Do NOT close the Terminal window that pops up when opening the App, it is required for it to run 

Getting Started 

1. In the white text box toward the top of the screen, enter the desired amino acid sequence. Be sure to type each amino acid using it corresponding letter/symbol as listed in the Legend. (For reference, you can click “Open Legend”, which is scrollable) 

2. Select the N- and C-terminus from the drop-down menu to the left and right of the text box. (N on the left and C on the right) 

3. Click “Generate Peptide” or press the ‘Enter’ key to populate the functions under all tabs. The chemical structure, mass, and net charge of the peptide should appear on the first page of the application. 

4. To enter a new amino acid sequence, clear the white text box at the top and type in the new sequence. Click “Generate Peptide”. 

Mass Spectrometry 

1. Click the “Mass Spectrometry” tab. 

2. On the left side of the tab, the mass spectrometry graph will appear and to the right a table of all the values produced at each peak. 

Reagent Calculations 

1. Click the “Conjugation” tab. 

2. Select whether you are going to use a dry or liquid reagent. The app will already default to dry.  

3. Input all necessary values in their corresponding boxes. Make sure the values are in the correct units (specified to the right of each input box). 

4. To calculate the ‘Reagent Mass’: 

	- Choose the dry option  

    - Then you must input numbers into fields 1 through 4 

    - If you want to calculate ‘Solvent Volume’, you must also input a number into the ‘Solvent Factor’ field 

5. To calculate the ‘Reagent Volume’: 

    - Choose the wet option 

    - Then you must input numbers into fields 1 through 5 

    - Again, if you want to calculate ‘Solvent Volume’, you must also input a number into the ‘Solvent Factor’ field 

Secondary Structure 

1. To view the predicted secondary structure, click the “Structure” tab. 

Adding Custom Amino Acid 

1. Click “Open Legend” at the top of the top of the application. 

2. At the bottom of the legend window, click “Add a new Amino Acid”. 

4. When inputting text into the fields remember: 

	- The Single Letter Code must be a single character that does not already exist in the database. 

		- Characters can be anything, including but not limited to letters, numbers, or symbols (!, ?, \, etc...) 

	- The Name or the SMILES of the Amino Acid cannot already exist in the database 

	- The Charge must be an integer value 

	- The Residue Mass must also be a floating-point number (which is a number that includes a decimal. 

5. If the inputs are not accepted, the program will report back an error when trying to add it.  

6. If the inputs are all valid, you will get a success message after clicking the button.  

7. The legend window must be closed and reopened to see the amino acids added. 
        """
        
        # Create a Text widget for the tutorial content
        tutorial_content = tk.Text(tutorial, wrap="word", font=("Arial", 12), bg="light grey", borderwidth=2, relief="solid")
        tutorial_content.insert("1.0", tutorial_text)

        # Define tags
        tutorial_content.tag_configure("list_number", font=("Arial", 14, "bold"))
        tutorial_content.tag_configure("heading_bold", font=("Arial", 12, "bold"))

        # Find all occurrences of the numbers followed by a period for list numbering
        import re
        for match in re.finditer(r"\d+\.", tutorial_text):
            start = "1.0 + {} chars".format(match.start())
            end = "1.0 + {} chars".format(match.end())
            tutorial_content.tag_add("list_number", start, end)

        # Apply bold tag to specified headings and sentences
        headings_and_sentences = [
            "NOTE: Do NOT close the Terminal window that pops up when opening the App, it is required for it to run",
            "Getting Started",
            "Mass Spectrometry",
            "Reagent Calculations",
            "Secondary Structure",
            "Adding Custom Amino Acid"
        ]
        for text in headings_and_sentences:
            start_idx = tutorial_content.search(text, "1.0", tk.END)
            if start_idx:
                end_idx = f"{start_idx} + {len(text)}c"
                tutorial_content.tag_add("heading_bold", start_idx, end_idx)

        tutorial_content.config(state="disabled")  # Make the text widget read-only
        tutorial_content.pack(padx=10, pady=10, fill="both", expand=True)

def open_legend():
    if not hasattr(open_legend, "is_open") or not open_legend.is_open:
        legend = Toplevel(app)
        legend.title("Amino Acid Legend")
        legend.minsize(300, 600)
        
        # Frame for Treeview and Scrollbar
        tree_frame = tk.Frame(legend)
        tree_frame.pack(expand=True, fill='both')
        
        # Create the treeview widget within the frame
        columns = ("Code", "Name")
        tree = ttk.Treeview(tree_frame, columns=columns, show="headings")
        tree.heading("Code", text="Input Code")
        tree.heading("Name", text="Name")
        
        # Adjust the columns' width to the content
        tree.column("Code", minwidth=50, width=100, anchor='center')
        tree.column("Name", minwidth=100, width=150, anchor='center')
        
        # Create a vertical scrollbar for the Treeview
        v_scroll = ttk.Scrollbar(tree_frame, orient="vertical", command=tree.yview)
        tree.configure(yscrollcommand=v_scroll.set)
        
        # Pack the treeview and scrollbar in the frame
        tree.pack(side=tk.LEFT, expand=True, fill='both')
        v_scroll.pack(side=tk.RIGHT, fill='y')

        # Read the CSV file using pandas and populate the treeview with the amino acids
        df = pd.read_csv("AminoAcidTable.csv")
        for index, row in df.iterrows():
            tree.insert("", tk.END, values=(row["Code"], row["Name"]))
        
        # Define the Add Amino Acid button with a command that includes itself, place it outside the frame
        add_AA_button = ttk.Button(legend, text="Add a new Amino Acid", command=lambda: open_add_AA_window(add_AA_button, legend))
        add_AA_button.pack(pady=10, side=tk.BOTTOM)
        
        legend.protocol("WM_DELETE_WINDOW", lambda: on_close_legend(legend))
        open_legend.is_open = True
        open_legend.add_AA_window = None  # Initially, there's no add_AA_window
        legend_button["state"] = "disabled"
        
def open_add_AA_window(button, legend):
    if not hasattr(legend, "add_AA_window_open") or not legend.add_AA_window_open:
        add_AA_window = Toplevel(app)
        add_AA_window.title("Add Amino Acids")
        add_AA_window.minsize(625, 350)

        button["state"] = "disabled"
        legend.add_AA_window_open = True  # Mark the add_AA_window as open
        legend.add_AA_window_ref = add_AA_window 

        # Use a main frame to help with centering content
        custom_aa_frame = ttk.Frame(add_AA_window)
        custom_aa_frame.place(relx=0.5, rely=0.5, anchor=tk.CENTER)

        # Heading Labels of this window
        add_AA_first_heading_label = tk.Label(custom_aa_frame, text="Add a custom amino acid here!", font=("Roboto", 18, "bold"))
        add_AA_first_heading_label.grid(row=0, column=0, columnspan=2, sticky=tk.N)
        add_AA_second_heading_label = tk.Label(custom_aa_frame, text="Input the characteristics of the Amino Acid below.", font=("Roboto", 12))
        add_AA_second_heading_label.grid(row=1, column=0, columnspan=2, sticky=tk.N, padx=15, pady=15)

        # Input labels and entry boxes
        custom_aa_input_fields = ["Single letter Code:", "Name of the Amino Acid:", "SMILES code:", "Charge:", "Residue Mass:"]
        entries = {}
        for i, field in enumerate(custom_aa_input_fields):
            label = ttk.Label(custom_aa_frame, text=f"{field}", font=('Helvetica', 12))
            label.grid(row=i + 2, column=0, sticky=tk.E, pady=5)
            entry = tk.Entry(custom_aa_frame)
            entry.grid(row=i + 2, column=1, sticky=tk.EW, padx=5, pady=5)
            entries[field.replace(":", "").replace(" ", "_").lower()] = entry  # Unique name assignment

        # Empty Label for success or failure message
        msg_label = tk.Label(custom_aa_frame, text="")
        msg_label.grid(row=7, column=0, columnspan=2, pady=2)

        # Create and place an "Add to Database" Button
        add_button = ttk.Button(custom_aa_frame, text="Add to Amino Acid Database", style="Modern.TButton", command=lambda: add_to_aaTable(entries, msg_label))
        add_button.grid(row=8, column=0, columnspan=2, sticky=tk.EW, padx=5, pady=5)

        add_AA_window.protocol("WM_DELETE_WINDOW", lambda: on_close_add_AA(add_AA_window, button, legend))


def add_to_aaTable(entries, msg_label, csv_path="AminoAcidTable.csv"):
    # Read the existing amino acid data
    try:
        aa_data = pd.read_csv(csv_path)
    except FileNotFoundError:
        msg_label.config(text="Error: Amino Acid database not found.", foreground="red")
        return
    
    # Extract the input data from entries
    data = {name: entry.get().strip() for name, entry in entries.items()}

    ## Input Validation
    # Validation 1: Check all inputs are filled and not just spaces
    if any(not value or value.isspace() for value in data.values()):
        msg_label.config(text="Error: All fields must be filled, cannot be blank or just spaces.", foreground="red")
        return

    # Validation 2: Check for spaces in all input boxes
    if any(' ' in value for value in data.values()):
        msg_label.config(text="Error: Inputs cannot contain spaces.", foreground="red")
        return
    
    # Validation 3: Check the Single Letter Code is a single character
    if len(data.get("single_letter_code", "")) != 1:
        msg_label.config(text="Error: Single letter code must be a single character.", foreground="red")
        return

    # Validation 4: Check types for charge and residue mass
    try:
        data['charge'] = int(data.get('charge', 0))
        data['residue_mass'] = float(data.get('residue_mass', 0.0))
    except ValueError:
        msg_label.config(text="Error: Charge must be an integer and residue mass must be a float.", foreground="red")
        return
    
    # Validation 5: Check if the single_letter_code already exists in the CSV
    if data["single_letter_code"] in aa_data["Code"].values:
        msg_label.config(text="Error: Amino Acid with this single letter code already exists.", foreground="red")
        return
    
    # Validation 6: Check if the Name of the Amino Acid already exists in the CSV
    if data["name_of_the_amino_acid"].upper() in aa_data["Name"].str.upper().values:
        msg_label.config(text="Error: Amino Acid with this name already exists.", foreground="red")
        return
    
    # Validation 7: Check if the SMILES of the Amino Acid already exists in the CSV
    if data["smiles_code"] in aa_data["SMILES"].values:
        msg_label.config(text="Error: Amino Acid with this SMIlES already exists.", foreground="red")
        return

    # All validation passed, then it appends the inputs into the csv
    new_row = pd.DataFrame([{
        "Code": data["single_letter_code"].upper(),
        "Name": data["name_of_the_amino_acid"],
        "SMILES": data["smiles_code"],
        "Charge": int(data["charge"]),
        "Residue Mass": float(data["residue_mass"])
    }])
    
    # Append new row to the DataFrame
    aa_data = pd.concat([aa_data, new_row], ignore_index=True)
    aa_data.to_csv(csv_path, index=False)
    msg_label.config(text="Amino Acid Added Successfully!", foreground="green")


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
app.geometry("1100x725")
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
error_msg = tk.Label(app, background="white", foreground= 'red', text="", font=("Roboto", 10))
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

# Can insert a group pic into the widget
# image_path = "Images\dog temp image.gif"  # Replace this with the path to your image file
# group_photo = PhotoImage(file=image_path)
# credits_info.image_create("end", image=group_photo)  # This inserts the image at the end of the text

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
