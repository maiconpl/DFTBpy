DFTBpy
======

DFTBpy is an open source project and the result of a long spare time work and effort of grokking and implementing the DFTB/SCC-DFTB method and pushing the limits of Python for speeding it up. Our focus is in achieving a balance between Software Engineering, easy maintenance/extension and a roughly good performance.

The aims of the project are:

1) To have an easier to modify and extend computational chemistry code based on Density Functional Theory (DFT): SCC-DFTB. The program is developed in Python with C extensions, for performance reasons.

2) To motivate undergraduate and graduate students as well as researchers to participate of the software development cycle.

3) To allow the students to understand, due to the nice features of Python, how a computational chemistry/physics
   program works. We believe the collaborators can participate of the project, extend it and mainly understand other computational chemistry methods and features.

4) To motivate people to learn and use Python, for scientific programming, as well as the Object Oriented paradigm.
   With Python, you can improve your research in computational chemistry/physics through: understanding computational methods (softwares written in Python are easier to understand), managing, automatically, several data got from simulations...You can develop softwares in Python for this.

5) A lot of computational chemistry/physics softwares, nowadays, have been developed in Python.
   Python has, though, achieved the scientific community in a way, as never seen before, since the advent of FORTRAN.

6) The collaborators and users of DFTBpy is welcome to ask any question about the
   software, the project and the DFTB/SCC-DFTB methodology by itself.

7) To perform research in Hight Performance Computing (HPC) using an interpreted and dynamic programming language.

8) You can contribute to this project in several ways:

   8.1) Giving us feedback about bugs and possible restrictions of the program.
   
   8.2) Testing it.
   
   8.3) Helping programming it.
   
   8.4) All critics are welcome.
   
   8.5) ...
   
9) A lot features, from programming to SCC-DFTB extensions, are going to be implemented:

   9.1) Optimize even more the code.
   
   9.2) Yet for speed purposes, port part of the demanding code to C:
        the Slater-Koster transformation and the creation of the H and S matrices.
        
   9.3) Explore parallelism in Python through the possible programming models: threads, MPI and, in the future, PyCUDA.
   
   9.4) Implement other geometry optimization algorithms: Gradient Conjugate (CG), for instance.
   
   9.5) Periodic boundary condition (PBC).
   
   9.6) DFTB/SCC-DFTB dynamic.

   9.7) Push the limits of Python.

# CREDITS:

The leader of the project (coder): Dr. Maicon Pierre Lourenco (MPL).

Collaborator: Dr. Mauricio Chagas da Silva.

Research group: GPQIT. zeus.qui.ufmg.br/~duarteh

This code can be downloaded, distributed and modified for FREE under the GPL (GNU General Public License) license terms.

# DISCLAIMER:

With this package distribution, only the SLAKO Mg-Mg.skf is distributed, with the aim of testing the code installation.
It was developed in the GPQIT group. To have access to more SLAKOs, please, see the www.dftb.org.

# ACKNOWLEDGMENT:

- Professor Helio Anderson Duarte

- GPQIT research group in Brazil:
  zeus.qui.ufmg.br/~duarteh

- Brazilian financial support:
  FAPEMIG, CAPES, CNPq, INCT Acqua.

- Collaborator:

  Dr. Mauricio Chagas da Silva

