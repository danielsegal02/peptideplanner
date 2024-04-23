import tkinter as tk
from tkinter import PhotoImage, ttk, Toplevel 
from SecondaryStructure.SecStructConnector import generate_sec_struct, retrieve_input
import os 

def display_text():
    # Get the text from the entry box and print it
    input_text = entry.get()
    generate_sec_struct(input_text)
    print(input_text)

# Create the main window
root = tk.Tk()
root.geometry("500x300")
root.title("Secondary Structure")

# Create an Entry widget (to input text)
entry = tk.Entry(root, width=40)
entry.pack(pady=20)  # Add some vertical padding

# Create a Button widget
button = tk.Button(root, text="Generate Structure", width = 15, command=display_text)
button.pack(pady=1)  # Add some vertical padding


tutorial_text = """
Secondary Structure Generation
Note: There is no need to re-enter the amino acid sequence, it will be automatically sent to the application.
1. Press the Generate Structure button, a tab in your browser will open, this may take a few seconds.
2. A code will appear on the screen, do not edit the code. 
3. At the top of the bar, under the line starting with file there will be an image with a save logo. Follow this line to the right until you reach the double arrow icon. Click this once.
4. The code will execute and you can scroll down to Cell 4, the cells are labeld on the left side with numbers. 
5. The image will generate there, it is 3D and it is rotatable, find the right angle and take a screenshot with the snipping tool built into your PC.

                        """
def display_tut():
    new_window = tk.Toplevel(root)
    text_label = tk.Label(new_window, text= tutorial_text, justify = "left")
    text_label.pack()

# Create a button
button2 = tk.Button(root, text="Tutorial", width = 15, command=display_tut)
button2.pack(pady = 1)


root.mainloop()
