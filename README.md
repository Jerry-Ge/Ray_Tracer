# Basic Ray Tracer
> **Author: Yuzhou Ge** -University of Southern California

### Operating System: 
	Max OS


### Main Features:
                               
-------------------------------------    
	1) Ray tracing triangles                  

	2) Ray tracing sphere                     

	3) Triangle Phong Shading                

	4) Sphere Phong Shading                   

	5) Shadows rays                           

	6) Still images                           
   
	7) Extra
	
	Recursive reflection
	I did the recursive reflection in function "traceRay". In order to use this function: when running the program in command line, add one more argument"-OR" in the end. !!!!!!!!!!!
	E.g: "./assign3 test.scene.txt test.jpg -OR"

	Good antialiasing
	I also have a good antialiasing by using the technqiue of super sampling. 
	You can check from the still images, the boundaries are smooth.

### Still Images:
	001: from screenfile.txt, no antialiasing applied, no recursive reflection
	002: from screenfile.txt, with antialiasing applied, no recursive reflection
	003: from test2.scene.txt, with antialiasing applied, no recursive reflection
	004: from test2.scene.txt, with antialiasing applied, with recursive reflection
	005: from spheres.scene.txt, with antialiasing applied, no recursive reflection
	006: from spheres.scene.txt, with antialiasing applied, with recursive reflection
	007: from table.scene.txt, with antialiasing applied, no recursive reflection
	008: from SIGGRAPH.scene.txt, with antialiasing applied, no recursive reflection