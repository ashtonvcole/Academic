import tkinter as tk
from tkinter import ttk
from tkinter import filedialog
from tkinter import messagebox
from tkinter import scrolledtext as st





class StartupWindow(tk.Frame):
    def __init__(self, parent, *args, **kwargs):
        super().__init__(parent, *args, **kwargs)
        self.parent = parent
        self.pack()

        # Variables
        self.nd = tk.StringVar()
        self.ele = tk.StringVar()
        self.disp = tk.StringVar()
        self.f = tk.StringVar()
        self.out = tk.StringVar()
        self.solve_ured = tk.IntVar(value = 1)
        self.solve_F_l = tk.IntVar()
        self.solve_epsilon = tk.IntVar()
        self.solve_sigma = tk.IntVar()
        self.solve_F_int = tk.IntVar()
        self.solve_yield = tk.IntVar()

        # Set up file selection window

        self.top_frame = tk.Frame()
        self.info_label = ttk.Label(self.top_frame, text = 'Select Truss Case')
        self.dir_button = tk.Button(self.top_frame, text = 'Select Directory', command = self.set_chad)

        self.files_frame = tk.Frame()
        self.nd_entry = tk.Entry(self.files_frame, state = 'disabled', width = 100, textvariable = self.nd)
        self.nd_button = tk.Button(self.files_frame, text = 'Select Nodes', command = self.set_nd)
        self.ele_entry = tk.Entry(self.files_frame, state = 'disabled', width = 100, textvariable = self.ele)
        self.ele_button = tk.Button(self.files_frame, text = 'Select Elements', command = self.set_ele)
        self.disp_entry = tk.Entry(self.files_frame, state = 'disabled', width = 100, textvariable = self.disp)
        self.disp_button = tk.Button(self.files_frame, text = 'Select Displacements', command = self.set_disp)
        self.f_entry = tk.Entry(self.files_frame, state = 'disabled', width = 100, textvariable = self.f)
        self.f_button = tk.Button(self.files_frame, text = 'Select Forces', command = self.set_f)
        self.out_entry = tk.Entry(self.files_frame, state = 'disabled', width = 100, textvariable = self.out)
        self.out_button = tk.Button(self.files_frame, text = 'Select Output Directory', command = self.set_out)

        self.options_frame = tk.Frame()
        
        self.solve_frame = tk.Frame(self.options_frame)
        self.solve_frame_label = ttk.Label(self.solve_frame, text = 'Solve for...')
        self.solve_ured_check = tk.Checkbutton(self.solve_frame, text = 'Reactionary Displacements', variable = self.solve_ured, state = 'disabled', onvalue = 1, offvalue = 0, command = self.toggle_ured)
        self.solve_F_l_check = tk.Checkbutton(self.solve_frame, text = 'Reactionary Forces', variable = self.solve_F_l, state = 'normal', onvalue = 1, offvalue = 0, command = self.toggle_F_l)
        self.solve_epsilon_check = tk.Checkbutton(self.solve_frame, text = 'Internal Strains', variable = self.solve_epsilon, state = 'normal', onvalue = 1, offvalue = 0, command = self.toggle_epsilon)
        self.solve_sigma_check = tk.Checkbutton(self.solve_frame, text = 'Internal Stresses (Average)', variable = self.solve_sigma, state = 'normal', onvalue = 1, offvalue = 0, command = self.toggle_sigma)
        self.solve_F_int_check = tk.Checkbutton(self.solve_frame, text = 'Internal Forces', variable = self.solve_F_int, state = 'normal', onvalue = 1, offvalue = 0, command = self.toggle_F_int)
        
        self.critical_frame = tk.Frame(self.options_frame)
        self.solve_yield_check = tk.Checkbutton(self.critical_frame, text = 'Identify Critical Yield, Crushing, and Buckling Elements', variable = self.solve_yield, state = 'normal', onvalue = 1, offvalue = 0, command = self.toggle_solve_yield)
        self.critical_frame_label = ttk.Label(self.critical_frame, text = 'This requires that the yield strength, \ncrushing strength, and second moment \nof area be provided for each element. \nAdditionally, it provides a comparison \nof the actual and critical values, as well \nas the minimum applied force \ndistribution causing failure as a \nmultiple of the supplied distribution.')
        
        self.execution_frame = tk.Frame()
        self.run_button = tk.Button(self.execution_frame, text = 'Solve', command=self.run, state = 'disabled')
        self.quit_button = tk.Button(self.execution_frame, text = 'Cancel', command = quit)

        # Set up file selection elements

        parent.title('Truss Solver ND - Startup')

        self.top_frame.pack()
        self.info_label.pack(pady = 10)
        self.dir_button.pack(pady = 10)

        self.files_frame.pack()
        self.files_frame.columnconfigure(0, weight = 3)
        self.files_frame.columnconfigure(1, weight = 1)
        self.nd_entry.grid(row = 0, column = 0, padx = 10, pady = 5)
        self.nd_button.grid(row = 0, column = 1, padx = 10, pady = 5, sticky = 'W')
        self.ele_entry.grid(row = 1, column = 0, padx = 10, pady = 5)
        self.ele_button.grid(row = 1, column = 1, padx = 10, pady = 5, sticky = 'W')
        self.disp_entry.grid(row = 2, column = 0, padx = 10, pady = 5)
        self.disp_button.grid(row = 2, column = 1, padx = 10, pady = 5, sticky = 'W')
        self.f_entry.grid(row = 3, column = 0, padx = 10, pady = 5)
        self.f_button.grid(row = 3, column = 1, padx = 10, pady = 5, sticky = 'W')
        self.out_entry.grid(row = 4, column = 0, padx = 10, pady = 5)
        self.out_button.grid(row = 4, column = 1, padx = 10, pady = 5, sticky = 'W')

        self.options_frame.pack()

        self.solve_frame.grid(row = 0, column = 0, padx = 10, pady = 5)
        self.solve_frame_label.pack()
        self.solve_ured_check.pack()
        self.solve_F_l_check.pack()
        self.solve_epsilon_check.pack()
        self.solve_sigma_check.pack()
        self.solve_F_int_check.pack()

        self.critical_frame.grid(row = 0, column = 1, padx = 10, pady = 5)
        self.solve_yield_check.grid(row = 0, column = 0)
        self.critical_frame_label.grid(row = 1, column = 0)

        self.execution_frame.pack(side = tk.RIGHT)
        self.quit_button.grid(row = 0, column = 3, padx = 5, pady = 5)
        self.run_button.grid(row = 0, column = 4, padx = 5, pady = 5)

        self.determine_if_solveable()
    
    def determine_if_solveable(self):
        if self.nd.get() != '' \
            and self.ele.get() != '' \
                and self.disp.get() != '' \
                    and self.f.get() != '' \
                        and self.out.get() != '':
            self.run_button['state'] = 'normal'
        else:
            self.run_button['state'] = 'disabled'
    
    def set_chad(self):
        folder = filedialog.askdirectory()
        if folder == '':
            return
        self.nd.set(f'{folder}/nodes')
        self.ele.set(f'{folder}/elements')
        self.disp.set(f'{folder}/displacements')
        self.f.set(f'{folder}/forces')
        self.out.set(folder)
        self.determine_if_solveable()
    
    def set_nd(self):
        path = filedialog.askopenfilename()
        if path == '':
            return
        self.nd.set(path)
        self.determine_if_solveable()

    def set_ele(self):
        path = filedialog.askopenfilename()
        if path == '':
            return
        self.ele.set(path)
        self.determine_if_solveable()

    def set_disp(self):
        path = filedialog.askopenfilename()
        if path == '':
            return
        self.disp.set(path)
        self.determine_if_solveable()
    
    def set_f(self):
        path = filedialog.askopenfilename()
        if path == '':
            return
        self.f.set(path)
        self.determine_if_solveable()
    
    def set_out(self):
        folder = filedialog.askdirectory()
        if folder == '':
            return
        self.out.set(folder)
        self.determine_if_solveable()

    def toggle_ured(self):
        pass

    def toggle_F_l(self):
        pass

    def toggle_epsilon(self):
        pass

    def toggle_sigma(self):
        pass

    def toggle_F_int(self):
        pass

    def toggle_solve_yield(self):
        pass

    def get_options(self) -> tuple:
        return (self.solve_ured.get() == 1, \
                self.solve_F_l.get() == 1, \
                self.solve_epsilon.get() == 1, \
                self.solve_sigma.get() == 1, \
                self.solve_F_int.get() == 1, \
                self.solve_yield.get() == 1)

    def run(self):
        # Validation first
        # No validation needed
        self.parent.destroy()





class ProgressWindow(tk.Frame):
    def __init__(self, parent, *args, **kwargs):
        super().__init__(parent, *args, **kwargs)
        self.parent = parent
        self.pack()

        # Variables
        self.progress = tk.IntVar()

        # Set up elements
        self.main_frame = tk.Frame()
        self.info_label = ttk.Label(self.main_frame, text = 'Solving...')
        self.bar = ttk.Progressbar(self.main_frame, orient = tk.HORIZONTAL, length = 300, variable = self.progress)
        self.logger = st.ScrolledText(self.main_frame, wrap = tk.WORD, width = 40, height = 5)

        # Set up window
        parent.title("Truss Solver ND - Progress")

        self.main_frame.pack()
        self.info_label.pack(pady = 10)
        self.bar.pack(pady = 10)
        self.logger.pack(padx = 20, pady = 10)
    
    def update_progress(self, progress: float, log_message: str = "") -> None:
        """Updates the progress bar and log.
        
        Args:
            progress: A float from 0.0 to 100.0 representing the percent to
                which the bar should be filled. A value of 100.0 makes the bar
                appear empty.
            log_message: A string added to the log. It should end with a
                newline character if the log needs a new line.
        
        Returns:
            None
        """
        self.progress.set(progress)
        self.logger.insert(tk.INSERT, log_message)
        self.logger.see(tk.END)
        self.parent.update()





if __name__ == '__main__':
    root = tk.Tk()
    sw = StartupWindow(root)
    sw.mainloop()