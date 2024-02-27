
To run two programs, you need to install two libraries (gmp-6.2.1, libpaillier-0.8), which are located in the library folder. After downloading and unzipping the files, you can install them as follows.

     cd gmp-6.2.1/
     ./configure
     make
     make check
     sudo make install
   

     cd libpaillier-0.8/
     ./configure 
     make
     make check
     sudo make install
