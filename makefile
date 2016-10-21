jam:
	cls
	gcc SPP_link_cell_2.3.c SPP_force.c SPP_sub.c -o jam -std=c99
	jam -p 1 < file_control.txt