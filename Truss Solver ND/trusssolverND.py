# Truss Solver ND





import numpy
import tkinter as tk
from tkinter import messagebox
import mygui





def line_list(the_file) -> list:
    """Reads a line of a file and returns space-delimited values a list.
    
    Args:
        the_file: An open file object.

    Returns:
        A list of the space-delimited entries on the most recent line.
    """
    return the_file.readline().strip("\n").split()





def outer_index_of(the_list: list, item, inner_index: int) -> int:
    """Returns the top index of an item in a list of lists.
    
    For a list of lists, this function will return the outer index of an item
    at a specified inner index within one of the lists. If nothing is found,
    -1 is returned.
    
    Args:
        the_list: A list of lists to be searched.
        item: The value sought.
        inner_index: The index to search at within each inner list.
    
    Returns:
        The index of the list in which the item is found, i.e. the outer index.
    """
    for outer_index in range(0, len(the_list)):
        if the_list[outer_index][inner_index] == item:
            return outer_index
    return -1





def length_of(element, nd, n_dim):
    """Gives the length of an element.
    
    Args:
        element: A list with entries representing an element.
        nd: A list of node lists.
        n_dim: The number of spatial dimensions the element resides in."""
    nd_nums = [outer_index_of(nd, element[1], 0), \
                outer_index_of(nd, element[2], 0)]
    L = 0
    for dim in range(1, n_dim + 1):
        L += (nd[nd_nums[1]][dim] - nd[nd_nums[0]][dim]) ** 2
    L = numpy.sqrt(L)
    return L





def get_K(nd_1: list, nd_2:list, E:float, A:float, n_dim: int) -> numpy.ndarray:
    """Gives a local stiffness matrix for a given element.

    This supports any N spatial dimensions greater than or equal to 1.
    
    Args:
        nd_1: A list representing the properties of the first node, in the
            form [nd_id, x, y, ...], where the number of coordinates,
            naturally, matches the number of dimensions considered for the
            truss structure.
        nd_2: A list representing the properties of the second node, in the
            same form.
        E: The Young's modulus of the element.
        A: The cross-sectional area of the element.
        n_dim: The number of dimensions in which the element resides.
    
    Returns:
        A (n_dim * 2)x(n_dim * 2) local stiffness matrix.
    """
    if n_dim > 0:
        # Find L
        L = 0
        for dim in range(1, n_dim + 1):
            L += (nd_2[dim] - nd_1[dim]) ** 2
        L = numpy.sqrt(L)
        # Find direction cosines
        c = []
        for dim in range(1, n_dim + 1):
            c.append((nd_2[dim] - nd_1[dim]) / L)
        # Compute matrix sub-block
        S = numpy.zeros((n_dim, n_dim))
        for dim1 in range(1, n_dim + 1):
            for dim2 in range(1, n_dim + 1):
                S[dim1 - 1][dim2 - 1] = E * A * c[dim1 - 1] * c[dim2 - 1] / L
        # Compute final matrix
        K = numpy.concatenate((S, numpy.negative(S)), axis = 0)
        K = numpy.concatenate((K, numpy.negative(K)), axis = 1)
        return K
    else:
        messagebox.showerror(title = 'Truss Solver ND - Execution Aborted', message = f'Error: unsupported dimension {n_dim}.')
        quit(1)





def solve_case(nodes_filename: str, elements_filename: str, displacements_filename: str, forces_filename: str, out_directory: str, options: list = (True, False, False, False, False, False), pw: mygui.ProgressWindow = None):

    # Open files
    try:
        nd_file = open(nodes_filename)
    except:
        messagebox.showerror(title = 'Truss Solver ND - Execution Aborted', message = f'Error: Unable to open nodes file \"{nodes_filename}\".')
        quit(1)
    try:
        ele_file = open(elements_filename)
    except:
        messagebox.showerror(title = 'Truss Solver ND - Execution Aborted', message = f'Error: Unable to open elements file \"{elements_filename}\".')
        quit(1)
    try:
        disp_file = open(displacements_filename)
    except:
        messagebox.showerror(title = 'Truss Solver ND - Execution Aborted', message = f'Error: Unable to open displacements file \"{displacements_filename}\".')
        quit(1)
    try:
        f_file = open(forces_filename)
    except:
        messagebox.showerror(title = 'Truss Solver ND - Execution Aborted', message = f'Error: Unable to open forces file \"{forces_filename}\".')
        quit(1)

    # Read in data
    # Please see PDF for conventions

    #####
    ##### CHECKPOINT 0/10 READING IN NODES
    #####
    if pw != None:
        pw.update_progress(0.0, 'Reading in nodes...\n')
    
    # Nodes
    nd = [] # A list of lists nd[nd_num] = [nd_id, x, y, ...]
    pos_of = [] # A list of lists, [nd_num][dim-1] = value is position in global matrix, initially
    line = line_list(nd_file)
    n_nd = int(line[0]) # Number of nodes
    n_dim = int(line[1]) # Number of dimensions
    if not n_dim > 0:
        messagebox.showerror(title = 'Truss Solver ND - Execution Aborted', message = f'Error: unsupported dimension {n_dim}.')
        quit(1)
    for nd_num in range(0, n_nd):
        line = line_list(nd_file)
        nd.append([])
        nd[nd_num].append(int(line[0])) # note_id
        pos_of.append([])
        for dim in range(1, n_dim + 1):
            nd[nd_num].append(float(line[dim])) # coordinate
            pos_of[nd_num].append(n_dim * nd_num + dim - 1)

    #####
    ##### CHECKPOINT 1/10 READING IN ELEMENTS
    #####
    if pw != None:
        pw.update_progress(10.0, 'Reading in elements...\n')
    
    # Elements
    ele = [] # A list of lists ele[ele_num] = [ele_id, nd_id_1, nd_id_2, E, A]
    line = line_list(ele_file)
    n_ele = int(line[0]) # Number of elements
    for ele_num in range(0, n_ele):
        line = line_list(ele_file)
        ele.append([])
        ele[ele_num].append(int(line[0])) # ele_id
        ele[ele_num].append(int(line[1])) # node_id_1
        ele[ele_num].append(int(line[2])) # node_id_2
        ele[ele_num].append(float(line[3])) # E
        ele[ele_num].append(float(line[4])) # A
        if options[5]:
            ele[ele_num].append(float(line[5])) # sigma_critical_y
            ele[ele_num].append(float(line[6])) # sigma_critical_c
            ele[ele_num].append(float(line[7])) # I

    #####
    ##### CHECKPOINT 2/10 READING IN FORCES
    #####
    if pw != None:
        pw.update_progress(20.0, 'Reading in forces...\n')
    
    # Forces
    f = [] # A list of lists f[f_num] = [node_id, dim, F]
    line = line_list(f_file)
    n_f = int(line[0])
    for f_num in range(0, n_f):
        line = line_list(f_file)
        f.append([])
        f[f_num].append(int(line[0])) # node_id
        f[f_num].append(int(line[1])) # dim
        f[f_num].append(float(line[2])) # F
    
    #####
    ##### CHECKPOINT 3/10 READING IN DISPLACEMENTS
    #####
    if pw != None:
        pw.update_progress(30.0, 'Reading in displacements...\n')

    # Displacements
    disp = [] # A list of lists disp[disp_num] = [node_id, dim, u]
    line = line_list(disp_file)
    n_disp = int(line[0])
    n_dof = n_dim * n_nd # Number of degrees of freedom of system
    for disp_num in range(0, n_disp):
        line = line_list(disp_file)
        disp.append([])
        disp[disp_num].append(int(line[0])) # node_id
        disp[disp_num].append(int(line[1])) # dim
        disp[disp_num].append(float(line[2])) # u
        the_nd_id = disp[disp_num][0]
        the_nd_num = outer_index_of(nd, the_nd_id, 0) # Find the list number of a given id
        the_dim = disp[disp_num][1]
        the_pos = pos_of[the_nd_num][the_dim - 1]
        for nd_num in range(0, n_nd):
            for dim in range(1, n_dim + 1):
                if pos_of[nd_num][dim - 1] > the_pos:
                    pos_of[nd_num][dim - 1] -= 1 # Move everyting after it forward
        pos_of[the_nd_num][the_dim - 1] = n_dim * n_nd - 1 # Move it to the end
        n_dof -= 1 # Update the added constraint

    # Close files
    nd_file.close()
    ele_file.close()
    disp_file.close()
    f_file.close()

    # Solve
    
    #####
    ##### CHECKPOINT 4 SETTING EQUATIONS
    #####
    if pw != None:
        pw.update_progress(40.0, 'Setting equations...\n')

    # [K_red](u_red) = (F_red), solve (u_red)
    # pos's 0 through n_dof - 1
    # [K_L](u) = (F_L), solve (F_L)
    # pos's n_dof through n_nd
    # L refers to lower, since visually we separate these constrained
    # equations and move them to the bottom of our matrix.

    K_red = numpy.zeros((n_dof, n_dof))
    F_red = numpy.zeros((n_dof, 1))
    K_l = numpy.zeros((n_dim * n_nd - n_dof, n_dim * n_nd))
    u_l = numpy.zeros((n_dim * n_nd - n_dof, 1))
    N_mat = numpy.zeros((n_ele, n_dim * n_nd)) # Strains for each node

    # Apply forces to F_red
    for force in f:
        the_pos = pos_of[outer_index_of(nd, force[0], 0)][force[1] - 1]
        F_red[the_pos] += force[2]

    # Arrange u_l
    for displacement in disp:
        the_pos = pos_of[outer_index_of(nd, displacement[0], 0)][displacement[1] - 1]
        u_l[the_pos - n_dof] = displacement[2]

    # Apply elements to K_red and subtract known displacements/stiffnesses from F_red
    # Also arrange N_mat
    for ele_num in range(0, n_ele):
        element = ele[ele_num]
        nd_nums = [outer_index_of(nd, element[1], 0), \
                   outer_index_of(nd, element[2], 0)] # Useful for the upcoming loop
        K_loc = get_K(nd[nd_nums[0]], \
            nd[nd_nums[1]], \
            element[3], \
            element[4], \
            n_dim)
        # Calculate L
        L = 0
        for dim in range(1, n_dim + 1):
            L += (nd[nd_nums[1]][dim] - nd[nd_nums[0]][dim]) ** 2
        L = numpy.sqrt(L)
        for loc_nd_num_a in range(0, 2):
            for dim_a in range(1, n_dim + 1):
                # Row/equation number
                pos_a = pos_of[nd_nums[loc_nd_num_a]][dim_a - 1] # Get global from local
                loc_pos_a = n_dim * loc_nd_num_a + dim_a - 1 # Nice parallel to global
                for loc_nd_num_b in range(0, 2):
                    for dim_b in range(1, n_dim + 1):
                        # Column/variable number
                        pos_b = pos_of[nd_nums[loc_nd_num_b]][dim_b - 1] # Get global from local
                        loc_pos_b = n_dim * loc_nd_num_b + dim_b - 1 # Nice parallel to global
                        # If pos_a (row) < n_dof, add K_loc to K_red or F_red
                        if pos_a < n_dof:
                            # If pos_b (col) < n_dof, add K_loc[loc_pos_a][loc_pos_b] to K_red[pos_a][pos_b]
                            if pos_b < n_dof:
                                K_red[pos_a][pos_b] += K_loc[loc_pos_a][loc_pos_b]
                            # Else, subtract K_loc[loc_pos_a][loc_pos_b] * u_l[pos_b - n_dof] from F_red[pos_a]
                            else:
                                F_red[pos_a] -= K_loc[loc_pos_a][loc_pos_b] * u_l[pos_b - n_dof]
                        # Else add K_loc[loc_pos_a][loc_pos_b] to K_L[pos_a - n_dof][pos_b]
                        else:
                            K_l[pos_a - n_dof][pos_b] += K_loc[loc_pos_a][loc_pos_b]
                # Freeload this loop to set value in N_mat
                # Negative if loc_nd_num_a == 0
                if loc_nd_num_a == 0:
                    N_mat[ele_num][pos_a] = -(nd[nd_nums[1]][dim_a] - nd[nd_nums[0]][dim_a]) / (L ** 2)
                # Positive if loc_nd_num_a == 1
                elif loc_nd_num_a == 1:
                    N_mat[ele_num][pos_a] = (nd[nd_nums[1]][dim_a] - nd[nd_nums[0]][dim_a]) / (L ** 2)
    
    #####
    ##### CHECKPOINT 5 SOLVING DISPLACEMENTS
    #####
    if pw != None:
        pw.update_progress(50.0, 'Solving displacements...\n')

    # Solve for u_red
    u_red = numpy.linalg.solve(K_red, F_red)
    # Append u_red and u_l into u
    u = numpy.concatenate((u_red, u_l), axis = 0)
    
    #####
    ##### CHECKPOINT 6 SOLVING EXTERNAL FORCES
    #####
    if options[1]:
        if pw != None:
            pw.update_progress(60.0, 'Solving external forces...\n')

        # Solve for F_l
        F_l = numpy.matmul(K_l, u)
    else:
        F_l = None

    #####
    ##### CHECKPOINT 7 SOLVING ELEMENT STRAINS
    #####
    # If you want stresses or forces or critical, you need strains
    if options[2] or options[3] or options[4] or options[5]:
        if pw != None:
            pw.update_progress(70.0, 'Solving element strains...\n')

        # Solve for epsilon
        epsilon = numpy.matmul(N_mat, u)
    else:
        epsilon = None

    #####
    ##### CHECKPOINT 8 SOLVING INTERNAL STRESSES
    #####
    # If you want critical (yield), you need stresses
    if options[3] or options[5]:
        if pw != None:
            pw.update_progress(80.0, 'Solving internal stresses...\n')

        # Solve for sigma
        # sigma[ele_num] = E * epsilon[ele_num]
        sigma = numpy.zeros((n_ele, 1))
        for ele_num in range(0, n_ele):
            sigma[ele_num] = ele[ele_num][3] * epsilon[ele_num]

    #####
    ##### CHECKPOINT 9 SOLVING INTERNAL FORCES
    #####
    # If you want critical (buckling), you need forces
    if options[4] or options[5]:
        if pw != None:
            pw.update_progress(90.0, 'Solving internal forces...\n')

        # Solve for F_int
        # F_int[ele_num] = E * A * epsilon[ele_num]
        F_int = numpy.zeros((n_ele, 1))
        for ele_num in range(0, n_ele):
            F_int[ele_num] = ele[ele_num][3] * ele[ele_num][4] * epsilon[ele_num]
    
    #####
    ##### CHECKPOINT 9.5 CRITICAL ELEMENTS
    #####
    if options[5]:
        if pw != None:
            pw.update_progress(95.0, 'Finding critical elements...\n')

        # F_int_cr = pi^2 * E * I / L^2
        # Since our system is linear...
        # Doubling load distribution means
        # Doubling internal force
        # Thus, load_mult_compression = F_int_critical / F_int
        # Likewise, load_mult_tension = sigma_critical / sigma
        
        ele_num_y = 0 # Element number of critical yield element
        min_ele_sigma_y_critical = 0 # Critical
        min_ele_sigma_y = 0 # Observed
        min_load_mult_y = float('inf') # What factor of P?
        ele_num_c = 0 # Element number of critical crushing element
        min_ele_sigma_c_critical = 0 # Critical
        min_ele_sigma_c = 0 # Observed
        min_load_mult_c = float('inf') # What factor of P?
        ele_num_b = 0 # Element number of critical buckling element
        min_ele_F_int_critical = 0 # Critical
        min_ele_F_int = 0 # Observed
        min_load_mult_b = float('inf') # What factor of P?
        # Go through each element
        for ele_num in range(0, n_ele):
            # If tensile, calculate multiple for yield
            if sigma[ele_num] > 0:
                ele_sigma_y_critical = ele[ele_num][5]
                ele_sigma_y = sigma[ele_num]
                load_mult_y = ele_sigma_y_critical / ele_sigma_y
                if load_mult_y < min_load_mult_y:
                    ele_num_y = ele_num
                    min_ele_sigma_y_critical = ele_sigma_y_critical
                    min_ele_sigma_y = ele_sigma_y
                    min_load_mult_y = load_mult_y
            # If compressive, calculate multiples for crushing and buckling
            elif sigma[ele_num] < 0:
                # Crushing
                ele_sigma_c_critical = ele[ele_num][6]
                ele_sigma_c = sigma[ele_num]
                load_mult_c = ele_sigma_c_critical / ele_sigma_c
                if load_mult_c < min_load_mult_c:
                    ele_num_c = ele_num
                    min_ele_sigma_c_critical = ele_sigma_c_critical
                    min_ele_sigma_c = ele_sigma_c
                    min_load_mult_c = load_mult_c
                # Buckling
                L = length_of(ele[ele_num], nd, n_dim)
                ele_F_int_critical = - ((numpy.pi ** 2 * ele[ele_num][3] * ele[ele_num][7]) / (L ** 2))
                ele_F_int = F_int[ele_num]
                load_mult_b = ele_F_int_critical / ele_F_int
                if load_mult_b < min_load_mult_b:
                    ele_num_b = ele_num
                    min_ele_F_int_critical = ele_F_int_critical
                    min_ele_F_int = ele_F_int
                    min_load_mult_b = load_mult_b


    #####
    ##### CHECKPOINT 10 SAVING RESULTS TO FILE
    #####
    if pw != None:
        pw.update_progress(99.9, 'Saving results to file...\n')
    
    # Here's everything we found
    # Reactionary Displacements
    # Internal Strains
    # Internal Stresses
    # Internal Forces
    # Reactionary Forces

    # Reactionary Displacements and Forces
    u_red_file = open(f'{out_directory}/reactionary_displacements', 'w')
    u_red_file.write(f'{len(u_red)}\n')
    # This is interspaced to try to save some time
    if options[1]:
        F_l_file = open(f'{out_directory}/reactionary_forces', 'w')
        F_l_file.write(f'{len(F_l)}\n')
    for nd_num in range(0, n_nd):
        for dim in range(1, n_dim + 1):
            the_pos = pos_of[nd_num][dim - 1]
            # If not already constrained, write displacement to file
            if the_pos < n_dof:
                u_red_file.write(f'{nd[nd_num][0]} {dim} {float(u_red[the_pos])}\n')
            # Else write force to file
            else:
                if options[1]:
                    F_l_file.write(f'{nd[nd_num][0]} {dim} {float(F_l[the_pos - n_dof])}\n')
    u_red_file.close()
    if options[1]:
        F_l_file.close()

    # Internal strains, stresses, forces, and critical
    if options[2]:
        epsilon_file = open(f'{out_directory}/internal_strains', 'w')
        epsilon_file.write(f'{n_ele}\n')
        for ele_num in range(0, n_ele):
            epsilon_file.write(f'{ele[ele_num][0]} {float(epsilon[ele_num])}\n')
        epsilon_file.close()
    if options[3]:
        sigma_file = open(f'{out_directory}/internal_stresses', 'w')
        sigma_file.write(f'{n_ele}\n')
        for ele_num in range(0, n_ele):
            sigma_file.write(f'{ele[ele_num][0]} {float(sigma[ele_num])}\n')
        sigma_file.close()
    if options[4]:
        F_int_file = open(f'{out_directory}/internal_forces', 'w')
        F_int_file.write(f'{n_ele}\n')
        for ele_num in range(0, n_ele):
            F_int_file.write(f'{ele[ele_num][0]} {float(F_int[ele_num])}\n')
        F_int_file.close()
    if options[5]:
        critical_file = open(f'{out_directory}/critical', 'w')
        critical_file.write('Failure Type - Element ID - Critical Value - ObservedValue - Factor of Load Causing Failure\n')
        critical_file.write(f'YieldStress {ele[ele_num_y][0]} {float(min_ele_sigma_y_critical)} {float(min_ele_sigma_y)} {float(min_load_mult_y)}\n')
        critical_file.write(f'CrushingStress {ele[ele_num_c][0]} {float(min_ele_sigma_c_critical)} {float(min_ele_sigma_c)} {float(min_load_mult_c)}\n')
        critical_file.write(f'BucklingForce {ele[ele_num_b][0]} {float(min_ele_F_int_critical)} {float(min_ele_F_int)} {float(min_load_mult_b)}\n')
        critical_file.close()

    #####
    ##### PROCESS COMPLETE
    #####
    if pw != None:
        pw.update_progress(100.0, 'Process complete\n')





def main():
    swroot = tk.Tk()
    sw = mygui.StartupWindow(swroot)
    sw.mainloop()

    pwroot = tk.Tk()
    pw = mygui.ProgressWindow(pwroot)
    solve_case(sw.nd.get(), sw.ele.get(), sw.disp.get(), sw.f.get(), sw.out.get(), sw.get_options(), pw)
    pw.mainloop()





if __name__ == "__main__":
    main()