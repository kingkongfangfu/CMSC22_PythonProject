import tkinter as tk
import structureFunctions as sf
import structure2Genalex as s2g
import genalexFunctions as gf

import sys

from tkinter import filedialog
from tkinter import simpledialog, messagebox

# Terminal display class
class TerminalDisplay(tk.Frame):
    def __init__(self, master=None, **kwargs):
        super().__init__(master, **kwargs)
        self.text_area = tk.Text(self, wrap="word", height=15, bg="black", fg="white", insertbackground="white")
        self.text_area.pack(fill="both", expand=True)
        self.text_area.config(state="disabled")

    def print_text(self, text):
        self.text_area.config(state="normal")
        self.text_area.insert(tk.END, text)
        self.text_area.see(tk.END)
        self.text_area.config(state="disabled")

    def clear_text(self):
        self.text_area.config(state="normal")
        self.text_area.delete("1.0", tk.END)
        self.text_area.config(state="disabled")

# Redirector class to capture print and errors
class ConsoleRedirector:
    def __init__(self, terminal_display):
        self.terminal_display = terminal_display

    def write(self, message):
        if message.strip() != "":
            self.terminal_display.print_text(message)

    def flush(self):
        pass

# Create main window
window = tk.Tk()
window.title("StruK2Stat")
window.geometry("700x700")
window.configure(bg="lightblue")

# Title
title_label = tk.Label(window, text="StruK2Stat", font=("Helvetica", 24), bg="lightblue")
title_label.pack(pady=20)

# Instruction
label = tk.Label(window, text="Select a file to open", bg="lightblue")
label.pack(pady=10)

# Buttons
tk.Button(window, text="Open File", command=sf.openStructureFile).pack(pady=5)
tk.Button(window, text="Run Structure", command=sf.runStructure).pack(pady=5)
tk.Button(window, text="Print Structure File", command=lambda: print(sf.structureFile)).pack(pady=5)
tk.Button(window, text="Edit STRUCTURE Main Parameters", command=sf.editMainParameters).pack(pady=5)
tk.Button(window, text="Edit STRUCTURE Extra Parameters", command=sf.editExtraParameters).pack(pady=5)
tk.Button(window, text="Convert to GenAlEx", command=s2g.struc2Genalex).pack(pady=5)
tk.Button(window, text="Compute Population Statistics", command=gf.populationStatistics).pack(pady=5)

# Terminal output display
terminal_display = TerminalDisplay(window)
terminal_display.pack(fill="both", expand=True, padx=10, pady=10)

# Redirect stdout and stderr to the terminal display
sys.stdout = ConsoleRedirector(terminal_display)
sys.stderr = ConsoleRedirector(terminal_display)

# Initial message
print("Welcome to StruK2Stat terminal.")

# Start the GUI event loop
window.mainloop()
# Reset stdout and stderr to default

