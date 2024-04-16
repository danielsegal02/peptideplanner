import tkinter as tk
from tkinter import PhotoImage, ttk
from test_args import generate_sec_struct, retrieve_input
import os 

figure_path = "figure_b.png"

def display_text():
    # Get the text from the entry box and print it
    input_text = entry.get()
    generate_sec_struct(input_text)
    chem_struct_photo = PhotoImage(file="figure_b.png")
    chem_struct_image_label.config(image=chem_struct_photo)
    chem_struct_image_label.photo = chem_struct_photo
    print(input_text)

# Create the main window
root = tk.Tk()
root.geometry("1100x725")
root.title("Simple Tkinter App")

# Create an Entry widget (to input text)
entry = tk.Entry(root, width=40)
entry.pack(pady=20)  # Add some vertical padding

# Create a Button widget
button = tk.Button(root, text="Read Input", command=display_text)
button.pack(pady=10)  # Add some vertical padding

# Create a Notebook to hold all the tabs
notebook = ttk.Notebook(root)
notebook.pack(fill='both',expand=True,padx=10,pady=10)

### Tab 1 for image and mass/charge info
tab1 = ttk.Frame(notebook)
notebook.add(tab1,text='Mass & Charge')

## Chemical structure image generation
# Frame to hold the image
chem_struct_image_frame = tk.Frame(tab1)
chem_struct_image_frame.configure(background="white")
chem_struct_image_frame.pack(fill='both',expand=True,padx=10,pady=10)

# Initialize the chem_struct_image_label with an empty image or placeholder
chem_struct_placeholder_image = tk.PhotoImage() # A blank PhotoImage object
chem_struct_image_label = tk.Label(chem_struct_image_frame,image=chem_struct_placeholder_image)
chem_struct_image_label.photo = chem_struct_placeholder_image  # Keep a reference
chem_struct_image_label.pack(padx=10, pady=10)



# Start the application
root.mainloop()
