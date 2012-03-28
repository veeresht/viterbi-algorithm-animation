#   Copyright 2012 Veeresh Taranalli <veeresht@gmail.com> 
#
#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program.  If not, see <http://www.gnu.org/licenses/>.

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
import matplotlib.patches as mpatches
import matplotlib.text as txtpatches
from commpy.channelcoding import Trellis
import matplotlib.animation as animation
import copy
import time

def generate_grid(trellis, trellis_length):
    """ Private method """

    grid = np.mgrid[0.12:0.22*trellis_length:(trellis_length+1)*(0+1j), 
                    0.1:0.1+trellis.number_states*0.1:trellis.number_states*(0+1j)].reshape(2, -1)
        
    return grid

def generate_states(trellis, trellis_length, grid, state_order, state_radius, font, state_metrics):
    """ Private method """
    state_patches = []
    #state_patches_array = np.empty([trellis_length, trellis.number_states])

    for time_index in xrange(trellis_length):
        
        patch_list = []
        for state_count in reversed(xrange(trellis.number_states)):
            state_pos = state_count + trellis.number_states*time_index
            state_patch = mpatches.Circle(grid[:,state_pos], state_radius, 
                    color="#003399", ec="#cccccc", alpha=1, label="1")
            patch_list.append(state_patch)
            #state_patches_array[time_index][state_count] = state_patch
            text_patch = txtpatches.Text(grid[0, state_pos], 
                                         grid[1, state_pos]-0.02, 
                                         str(state_order[state_count]), ha="center", 
                                         family=font, size = 20, color="#FFFFFF")
            patch_list.append(text_patch)
            metric_patch = txtpatches.Text(grid[0, state_pos],
                                           grid[1, state_pos] + 1.1*state_radius, 
                                           '(' + str(state_metrics[time_index][state_order[state_count]]) +')', 
                                           ha="center", family=font, size=14, color="#000000")
            patch_list.append(metric_patch) 
        
        state_patches.append(patch_list)
    
    return state_patches

def generate_edges(trellis, trellis_length, grid, state_order, state_radius, edge_colors, prev_states):
    """ Private method """
    edge_patches = []
 
    mapping_matrix = dict()

    for current_time_index in xrange(trellis_length-1):
        grid_subset = grid[:,trellis.number_states * current_time_index:]
        edge_list = []
        for state_count_2 in reversed(xrange(trellis.number_states)):
            input_count = 0
            prev_state = prev_states[current_time_index][trellis.number_states - state_count_2 - 1]
            for state_count_1 in xrange(trellis.number_states):
                dx = grid_subset[0, state_count_2+trellis.number_states] - grid_subset[0,state_count_1] - 2*state_radius
                dy = grid_subset[1, state_count_2+trellis.number_states] - grid_subset[1,state_count_1]
                if np.count_nonzero(trellis.next_state_table[state_order[state_count_1],:] == state_order[state_count_2]):
                    found_index = np.where(trellis.next_state_table[state_order[state_count_1],:] ==
                                            state_order[state_count_2])
                    
                    if prev_state == state_order[state_count_1]:
                        e_color = edge_colors[found_index[0][0]]
                    else:
                        e_color = "#DFDFDF"

                    edge_patch = mpatches.FancyArrow(grid_subset[0,state_count_1]+state_radius, 
                            grid_subset[1,state_count_1], dx, dy, width=0.005,
                            length_includes_head = True, color = e_color)
                    edge_list.append(edge_patch)

                    mapping_matrix[str(current_time_index) + 
                                   str(state_order[state_count_1]) + 
                                   str(state_order[state_count_2])] = edge_patch
                    input_count = input_count + 1
        
        edge_patches.append(edge_list)
    
    return [edge_patches, mapping_matrix]

def generate_labels(trellis, trellis_length, grid, state_order, state_radius, font, rx_codeword):
    """ Private method """

    for state_count in xrange(trellis.number_states):
        for input_count in xrange(trellis.number_inputs):
            edge_label = str(input_count) + "/" + str(
                    trellis.output_table[state_order[state_count], input_count])
            plt.text(grid[0, state_count]-1.5*state_radius, 
                     grid[1, state_count]+state_radius*(1-input_count-0.7), 
                     edge_label, ha="center", family=font, size=14)
            
    for current_time_index in xrange(trellis_length-1):
        plt.text(grid[0, 3] + current_time_index*0.2 + 0.1,
                grid[1, 3] + 0.1, str(rx_codeword[current_time_index][0]) + str(rx_codeword[current_time_index][1]),
                ha="center", family=font, size=14, color="#000000")

def to_1D_list(x):
    
    a = []
    for x_list in x:
        for element in x_list:
            a.append(element)
    return a


def viterbi_visualize(trellis, rx_codeword, previous_states, state_metrics, 
                      state_sequence, trellis_length=2, state_order=None, 
                      state_radius=0.04, edge_colors=None):
    """ Plot the Viterbi Algorithm animation. """

    if edge_colors is None:
        edge_colors = ["#9E1BE0", "#06D65D"]

    if state_order is None:
        state_order = range(trellis.number_states)
    
    font = "sans-serif"
    fig = plt.figure()
    ax = plt.axes([0,0,1,1])
    trellis_patches = []

    state_order.reverse()

    # Generate the elements of the Trellis
    trellis_grid = generate_grid(trellis, trellis_length)
    
    state_patches = generate_states(trellis, trellis_length, trellis_grid, 
                                    state_order, state_radius, font, 
                                    state_metrics)
    
    [edge_patches, mapping_matrix] = generate_edges(trellis, trellis_length, 
                                                    trellis_grid, state_order, 
                                                    state_radius, edge_colors, 
                                                    previous_states)
    
    generate_labels(trellis, trellis_length, trellis_grid, state_order, state_radius, font, rx_codeword)

    # Create frames for the animation
    frames = []
    
    for time_index in xrange(trellis_length-1):

        for spatch_num in xrange(0, trellis.number_states):

            frames.append(edge_patches[time_index]
                                      [(spatch_num)*trellis.number_inputs:
                                       (spatch_num+1)*trellis.number_inputs]) 
            frames.append(state_patches[time_index+1][(spatch_num)*3:(spatch_num+1)*3]) 
            
    # Add the trellis elements to the plot
    for state_list in state_patches:
        for patch in state_list:
            ax.add_artist(patch)

    for edge_list in edge_patches:
        for patch in edge_list:
            ax.add_artist(patch)
    
    
    # Highlight the best survivor path as the last step
    # of the Viterbi Algorithm
    state_sequence_patches = []
    for time_index in xrange(trellis_length-1):

        current_state = state_sequence[time_index]
        next_state = state_sequence[time_index+1]
        state_sequence_patches.append(mapping_matrix[str(time_index) + str(current_state) + str(next_state)])
    
    for edge_patch in reversed(state_sequence_patches):

        copy_patch = copy.copy(edge_patch)
        copy_patch.set_color("#000000")
        frames.append([copy_patch])
    
    ax.set_xticks([])
    ax.set_yticks([])
    
    ani = animation.ArtistAnimation(fig, frames, interval=1000, blit=True, repeat=False)
    #ani.save('viterbi_algorithm.mp4', fps=1, clear_temp=True)

    plt.show()


# =============================================================================
# Example showing the usage of viterbi_visualize
# =============================================================================
# Define the number of memory elements 
# per input in the convolutional encoder
memory = np.array([2])

# Define the generator matrix of 
# the convolutional encoder 
# Entries are in octal format
g_matrix = np.array([[05, 07]])

# Create a trellis representation 
# from the encoder generator matrix 
trellis = Trellis(memory, g_matrix)

# Specify the number of time steps 
# in the trellis diagram
trellis_length = 5

# Specify the order in which states 
# should appear on the trellis diagram
state_order = [0, 2, 1, 3]

# Specify the colors for 0, 1 inputs
# '#FF0000' --> Red   (edge corresponding to input 0)
# '#00FF00' --> Green (edge corresponding to input 1)
bit_colors = ['#FF0000', '#00FF00'] 

# Viterbi Algorithm outputs
rx_codeword = [[1, 1], [1, 1], [1, 1], [0, 0]] 
previous_states = [[0, 0, 10, 10], [0, 0, 2, 2], [1, 1, 3, 3], [0, 1, 3, 3]]
state_metrics = [[0, 'inf', 'inf', 'inf'], [2, 'inf', 0, 'inf'], [4, 1, 2, 1], [1, 2, 3, 2], [1, 3, 2, 3]]
state_sequence = [0, 2, 1, 0, 0]

# Plot the trellis diagram
viterbi_visualize(trellis, rx_codeword, previous_states, state_metrics, 
                  state_sequence, trellis_length, state_order, 
                  edge_colors = bit_colors)
