# NonRigidSimulation

1) Implementation of the FEPR 

https://www.cs.utah.edu/~ladislav/dinev18FEPR/dinev18FEPR.html

The FEPR functions are all in the Mesh.cpp file 

2) Build 

Go to directory 
mkdir build 
cd build 
cmake ..
cd .. 
cmake --build build 

To run : 
./FEPR

3) Toggle FEPR 

Go to the main.cpp file and comment the line in the update function + uncomment the other two. 

4) Mesh changes 

If you change the mesh make sure to change the high to the ground (m_ground value in mesh.cpp and m_O which stand for origin which is the center of rotation)



