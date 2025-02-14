import os.path
import tkinter as tk
from tkinter import ttk
from tkinter import font
from tkinter import filedialog

# Default values for data
data = {'avg': 1.0, 'std': 0.1, 'num': 1000, 'den': 0.25, 'olp': 0.0, 'dst': 'gamma', 'pbc': False, 'sar': False,
        'dir': './Data/user_data'}


class HelpGUI:
    def __init__(self, parent):
        # Create a top-level window
        self.top = tk.Toplevel(parent)
        self.top.title("Help")
        self.top.geometry("375x365")  # Width x Height

        # Create the fonts
        self.title_font = font.Font(family='Helvetica', size=15, weight='bold')

        # Dictionary to hold the help topics and their descriptions
        self.help_info = {
            "About": "Welcome to foam_gen Help!\n\nThis program can be used to generate statistical ensembles representing foams. ",
            "Average": "\"Average\" refers to the the average ball/bubble radius size that the sample distribution uses to sample from.",
            "CV": "\"CV\" or \"Coefficient of Variation\" is a general term that describes the width of the distribution. This is very similar to a ",
            "Number": "",
            "Density": "",
            "Overlap": "",
            "Distribution": "",
            "Periodic Boundary": "",
            "Standardized Atomic Radii": "",
            "Output Directory": ""
        }

        # Create the left and right gaps
        ttk.Label(self.top, text=" ").grid(row=0, column=0, padx=5)

        # Create the title
        ttk.Label(self.top, text="Foam Gen Help", font=self.title_font).grid(row=1, columnspan=2, column=1)

        # Create the label for the combobox
        self.comboboxlabel = ttk.Label(self.top, text="Help Topic:")
        self.comboboxlabel.grid(row=2, column=1, sticky='w', pady=10, padx=10)

        # Create a Combobox to select the help topic
        self.topic_var = tk.StringVar()
        self.combobox = ttk.Combobox(self.top, textvariable=self.topic_var, state="readonly", width=19)
        self.combobox['values'] = list(self.help_info.keys())
        self.combobox.grid(row=2, column=2, padx=10, pady=10, sticky='e')
        self.combobox.current(0)  # Default to first item in the list
        self.combobox.bind("<<ComboboxSelected>>", self.update_description)

        # Text widget or Label for displaying the help description
        self.description = tk.Text(self.top, height=15, width=40)
        self.description.grid(row=3, column=1, padx=10, pady=10, columnspan=2)
        self.description.insert('end', self.help_info[self.topic_var.get()])
        self.description.config(state='disabled')  # Make the text widget read-only

    def update_description(self, event):
        # Update the text widget with the selected topic's description
        self.description.config(state='normal')  # Enable text widget to modify
        self.description.delete('1.0', 'end')  # Clear current text
        self.description.insert('end', self.help_info[self.topic_var.get()])  # Insert new description
        self.description.config(state='disabled')  # Set back to read-only


def settings_gui():
    # Function to collect values and update the data dictionary
    global data

    # Browse the folder request function
    def browse_folder(browse_request="Choose Output Directory"):
        my_folder = filedialog.askdirectory(title=browse_request)
        dir_var.set(my_folder)
        dir_var_name.set(my_folder if len(my_folder) < 47 else '...' + my_folder[-44:])

    def help_gui():
        HelpGUI(root)

    def apply_values():
        global data
        data = {
            "avg": float(avg_var.get()),
            "std": float(std_var.get()),
            "num": int(num_var.get()),
            "den": float(den_var.get()),
            "olp": float(olp_var.get()),
            "dst": dst_var.get().lower(),
            "pbc": pbc_var.get(),
            "sar": sar_var.get(),
            "dir": dir_var.get() if os.path.exists(dir_var.get()) else None
        }
        root.destroy()

    # Function to handle cancel
    def cancel():
        root.destroy()

    # Main window
    root = tk.Tk()
    root.title("FoamGen")
    root.geometry("375x485")  # Set the window size
    # Styles
    background_color = 'wheat'
    foreground_color = 'black'
    # root.configure(bg=background_color)

    # Variables for storing input
    avg_var = tk.StringVar(value=str(data['avg']))
    std_var = tk.StringVar(value=str(data['std']))
    num_var = tk.StringVar(value=str(data['num']))
    den_var = tk.StringVar(value=str(data['den']))
    olp_var = tk.StringVar(value=str(data['olp']))
    dst_var = tk.StringVar(value=data['dst'].capitalize())
    pbc_var = tk.BooleanVar(value=data['pbc'])
    sar_var = tk.BooleanVar(value=data['sar'])
    dir_var = tk.StringVar(value=data['dir'])
    dir_var_name = tk.StringVar(value="./foam_gen" + data['dir'][1:])
    
    # Style
    my_style = ttk.Style(root)
    my_style.configure("Custom.TLabel", background=background_color, foreground=foreground_color)

    # Fonts
    title_font = font.Font(family="Cooper Black", size=25, weight='bold')
    setting_labels_font = font.Font(family='Serif', size=10, weight='bold')
    settings_font = font.Font(family='Serif', size=10)

    # Create the title
    test = ttk.Label(root, text='Foam Gen', font=title_font)
    test.grid(row=0, column=1, columnspan=2, pady=15)

    # Padding on the right
    ttk.Label(root, text=' ').grid(column=0, padx=10)

    # Setting up the grid and padding
    options = {'padx': 10, 'pady': 5}  # Common options for padding

    # Average Entry
    ttk.Label(root, text="Average", font=setting_labels_font).grid(row=1, column=1, sticky='w', **options)
    ttk.Entry(root, textvariable=avg_var).grid(row=1, column=2, **options)

    # Standard Deviation Entry
    ttk.Label(root, text="CV", font=setting_labels_font).grid(row=2, column=1, sticky='w', **options)
    ttk.Entry(root, textvariable=std_var).grid(row=2, column=2, **options)

    # Number Entry
    ttk.Label(root, text="Number", font=setting_labels_font).grid(row=3, column=1, sticky='w', **options)
    ttk.Entry(root, textvariable=num_var).grid(row=3, column=2, **options)

    # Density Entry/Slider
    ttk.Label(root, text="Density", font=setting_labels_font).grid(row=4, column=1, sticky='w', **options)
    ttk.Entry(root, textvariable=den_var).grid(row=4, column=2, **options)

    # Overlap Entry/Slider
    ttk.Label(root, text="Overlap", font=setting_labels_font).grid(row=5, column=1, sticky='w', **options)
    ttk.Entry(root, textvariable=olp_var).grid(row=5, column=2, **options)

    # Distribution Dropdown
    ttk.Label(root, text="Distribution", font=setting_labels_font).grid(row=6, column=1, sticky='w', **options)
    dst_options = ["Gamma", "Log-Normal", "Weibull", "Normal", "Half-Normal", "Physical 1 (Devries)",
                   "Physical 2 (Gal-Or)", "Physical 3 (Lemelich)"]
    dst_menu = ttk.Combobox(root, textvariable=dst_var, values=dst_options, width=14, font=settings_font)
    dst_menu.grid(row=6, column=2, **options)
    dst_menu.current(dst_options.index(data['dst'].capitalize()))

    # Periodic Boundary Condition Checkbox
    ttk.Label(root, text="Periodic Boundary", font=setting_labels_font).grid(row=7, column=1, sticky='w', **options)
    ttk.Checkbutton(root, variable=pbc_var).grid(row=7, column=2, **options)

    # Standardize Atomic Radii Checkbox
    ttk.Label(root, text="Standardize Atomic Radii", font=setting_labels_font).grid(row=8, column=1, sticky='w', **options)
    ttk.Checkbutton(root, variable=sar_var).grid(row=8, column=2, **options)

    # Browse output directory
    ttk.Label(root, text=" ").grid(row=9, pady=8)
    ttk.Label(root, text="Output Directory:", font=setting_labels_font).grid(row=10, column=1, sticky='w', padx=10)
    ttk.Button(root, text='Browse', command=browse_folder).grid(row=10, column=2, padx=10)
    ttk.Label(root, textvariable=dir_var_name, font=settings_font).grid(row=11, column=1, columnspan=2, padx=10, sticky='w')

    # Buttons
    ttk.Label(root, text=" ").grid(row=12, pady=8)
    ttk.Button(root, text="Help", command=help_gui).grid(row=13, column=1, sticky='w', pady=5)
    ttk.Button(root, text="Cancel", command=cancel).grid(row=13, column=1, sticky='e', pady=5)
    ttk.Button(root, text="Create Foam", command=apply_values).grid(row=13, column=2, sticky='e', **options)

    # Run the GUI
    root.mainloop()

    return data


if __name__ == '__main__':
    print(settings_gui())
