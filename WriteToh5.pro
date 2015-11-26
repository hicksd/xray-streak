PRO writeToh5,filename,struct
;this writes a structure to an h5 file
;
;First, create an IDL structure variable. It is essential that the structure not contain pointers, complex numbers, or object references, since these are not allowed in HDF5 files.
;struct = {name:"coyote", age:25, rascal:1B, salary:200000D, girlfriends:FltArr(546)}
;To write a structure into an HDF5 file, you need to create compound data type. As with all HDF files, 
;we first create a file identifier, since all communication with the file is done via the file identifier.
;filename = 'hdf5testfile.h5'
fileID = H5F_CREATE(filename)
;Next, we create a datatype object by supplying the structure as an argument. This will indicate that we 
;wish to create a compound data type, with the types determined from the fields of the structure.
datatypeID = H5T_IDL_CREATE(struct)
;The next step is to create a simple (or, in this case, not so simple) data space. As an argument, we 
;require the dimensions of the dataspace. In this case, with a single structure variable, the dimension is 1.
dataspaceID = H5S_CREATE_SIMPLE(1)
;Now we can create the dataset at a specified location, given as the first argument in the function below. 
;The second argument is the name of the dataset in the file, and the third and fourth arguments are the data 
;type identifier and data space identifier, respectively.
datasetID = H5D_CREATE(fileID, 'data', datatypeID, dataspaceID)
;Finally, we can write the structure into the file, and close the file.
H5D_WRITE, datasetID, struct
H5F_CLOSE, fileID

RETURN
END
