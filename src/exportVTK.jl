export exportMeshAndModelToVTK,CreateVTKtimeSeriesCollection

if hasJOcTree
function exportMeshAndModelToVTK(fnameNoExt::String,mesh::OcTreeMesh,data...)
  #writeOctreeAsVTKunstructured(fnameNoExt,mesh,data...)
  
  #Converts an octree mesh to an unstructured mesh and writes the result to an xml file
  #with compressed binary data storage.
  #Input:
  #  fnameNoExt: file name of output xml file without extension
  #  mesh: jInv octree mesh, can be FEM or FV mesh.
  #  data: Data to be attached to mesh. An arbitrary number of data fields can be attached
  #to the mesh, each as a separate argument. Each field must be written as a tuple in 
  #the format (dataArray,name), where dataArray is a dense array holding the data, it 
  #can be a vector or a multidimensional array. name is a string that gives the name
  #the data field will have in paraview.
  
  #Convert OcTree mesh to vtk unstructured grid format (get node coordinates and 
  #cell connectivity
  pts,cells = convertOctreeToVTKunstructured(mesh)
  
  #Generate xml file with mesh data
  vtk = vtk_grid(fnameNoExt,pts,cells)
  
  #Extract data from varargs and add data fields to xml file
  nnodes = size(pts,2)
  for idat in data
    if ~(typeof(idat[2])<:AbstractString)
        error("exportMeshAndModelToVTK: data must be input as tuple->(data::Array,description::String[,weights::Array])")
    end
    if (length(idat) == 2)
      addVTKdataOcTree(vtk,idat[1],idat[2],nnodes,mesh)
    else
      addVTKdataOcTree(vtk,idat[1],idat[2],nnodes,mesh,idat[3])
    end
  end
  
  #Write xml file to disk and close it
  vtk_save(vtk)
  println("Wrote mesh to vtk file with data fields:")
  for idat in data
     println(idat[2])
  end
end

#-----------------------------------------------

function convertOctreeToVTKunstructured(mesh::OcTreeMesh)
  S = mesh.S
  ncells = mesh.nc
  n1,n2,n3 = mesh.n
  ic,jc,kc,bsz = find3(S)
  inn,jnn,knn,Snnum = getNodalNumberingVTKlocal(S)
  nnodes = nnz(Snnum.SV)
  cells = Array(MeshCell,ncells)
  pts   = zeros(Float32,3,nnodes)
  const celltype = VTKCellTypes.VTK_HEXAHEDRON
  
  #Get connectivity
  #=  Nodes within a cell numbered as shown (Converted to vtk hexahdron numbering for
      input to WriteVTK functions
     /5----/7
    / |   /|
    1----3 |
    | /6-|-/8
    |/   |/
    2----4
  =#
  for cIdx = 1:ncells
    inds    = Array(Int32, 8)
    regNd   = (kc[cIdx]-1)*(n1+1)*(n2+1) + (jc[cIdx]-1)*(n1+1) + ic[cIdx]  #regNd is the
                                                    #underlying regular mesh node number
                                                    #of node 1 in current cell.
    hc      = S.SV.nzval[cIdx]
    inds[1] = Snnum.SV[regNd]
    inds[2] = Snnum.SV[regNd+hc]
    inds[3] = Snnum.SV[regNd+hc*(n1+1)]
    inds[4] = Snnum.SV[regNd+hc*(n1+1)+hc]
    inds[5] = Snnum.SV[regNd+hc*(n1+1)*(n2+1)]
    inds[6] = Snnum.SV[regNd+hc*(n1+1)*(n2+1)+hc]
    inds[7] = Snnum.SV[regNd+hc*(n1+1)*(n2+1) + hc*(n1+1)]
    inds[8] = Snnum.SV[regNd+hc*(n1+1)*(n2+1) + hc*(n1+1)+hc]
    newInd  = [inds[1];inds[2];inds[4];inds[3];inds[5];inds[6];inds[8];inds[7]]
    cells[cIdx] = MeshCell(celltype,newInd)
  end
  
  #Get node coordinates
  pts = getNodalGrid(mesh)'
  
  return pts, cells
end

#-----------------------------------------------

function getNodalNumberingVTKlocal(S::SparseArray3D)
  #Same as getNodalNumbering function in JOcTree/getNodalNumbering.jl except
  #that this function returns 3D indices of nodes.
  
  m1,m2,m3 = S.sz
  i,j,k,bsz = find3(S); bsz = round(Int64,bsz)

  nind = [i        j       k;
          i        j       k+bsz;
          i        j+bsz   k;
          i        j+bsz   k+bsz;
          i+bsz    j       k;
          i+bsz    j       k+bsz;
          i+bsz    j+bsz   k;
          i+bsz    j+bsz   k+bsz ]

  Ni = sparse3(nind[:,1],nind[:,2],nind[:,3],nind[:,3],[m1+1,m2+1,m3+1])
  i,j,k = find3(Ni)
  N = sparse3(i,j,k,1:length(i), [m1+1,m2+1,m3+1]);
  return i,j,k,N
end

# function getCellNumberingVTKlocal(S::SparseArray3D)
#   #Same as getCellNumbering function in JOcTree/getCellNumbering.jl except
#   #that this function returns 3D indices of cells.
#   m1,m2,m3  = S.sz
#   i,j,k     = find3(S)
#   return i,j,k,sparse3(i,j,k,1:length(i),[m1,m2,m3]);
# end

#-----------------------------------------------

function addVTKdataOcTree{T<:Real}(
                   vtk::WriteVTK.DatasetFile,data::Array{T},dataName::AbstractString,
                   nn::Int64,mesh::OcTreeMesh,weights=[])
  ndata = length(data)
  nc = mesh.nc
  nf = sum(mesh.nf)
  ne = sum(mesh.ne)
  @assert ndata in (nn,3*nn,nc,3*nc,nf,ne) "Length of data array does not match mesh"
  
  if ndata in (nn,3*nn)
    vtk_point_data(vtk,data,dataName)
  elseif ndata in (nc,3*nc)
    vtk_cell_data(vtk,data,dataName)
  elseif ndata == nf
    error("Face data not implemented")
  elseif ndata == ne
    error("Edge data not implemented for OcTree meshes")  
  end
  return
end

end #hasJOCTree if block
#--------------------------------------------------


function exportMeshAndModelToVTK(meshfile::String, mesh::TensorMesh3D,data...)

  x,y,z = getXYZvtk(mesh)
  
  #Initialise mesh file
  vtk = vtk_grid(meshfile,x,y,z)
  
  #Extract data from varargs and add data fields to xml file
  for idat in data
    if ~(typeof(idat[2])<:AbstractString)
        error("exportMeshAndModelToVTK: data must be input as tuple->(data::Array,description::String[,weights::Array])")
    end
    if (length(idat) == 2)
      addVTKdataTensor(vtk,idat[1],idat[2],mesh)
    else
      addVTKdataTensor(vtk,idat[1],idat[2],mesh,idat[3])
    end
  end
  
  outfiles = vtk_save(vtk)
  #println("Saved mesh and data to:", outfiles...)
  println("Wrote mesh to vtk file with data fields:")
  for idat in data
     println(idat[2])
  end
end

function addVTKdataTensor{T<:Real}(
                   vtk::WriteVTK.DatasetFile,data::Array{T},dataName::AbstractString,
                   mesh::TensorMesh3D,weights=[])
  nc = mesh.nc
  nn = prod(mesh.n + 1)
  nf = sum(mesh.nf)
  ne = sum(mesh.ne)
  ndata = length(data)
  @assert ndata in (nn,3*nn,nc,3*nc,nf,ne) "Length of data array does not match mesh"
  
  if ndata in (nn,3*nn)
    vtk_point_data(vtk,data,dataName)
  elseif ndata in (nc,3*nc)
    vtk_cell_data(vtk,data,dataName)
  elseif ndata == ne
    dcc = convertEdgeVecToCellCentreVec(data,mesh,weights)
    vtk_cell_data(vtk,dcc,dataName)
  elseif ndata == nf
    println("Face data not yet implemented")
  end
  return
end

function convertEdgeVecToCellCentreVec(data,mesh::TensorMesh3D,weights=[])
  if isempty(weights)
    weights = ones(mesh.nc)
  elseif (~isempty(weights) & (mesh.nc != length(weights)) )
    error("convertEdgeVecToCellCentreVec: data and weight lengths not equal")
  end
  wMat = spdiagm(weights)
  A1 = kron(av(mesh.n[3]),kron(av(mesh.n[2]),speye(mesh.n[1])))
  A2 = kron(av(mesh.n[3]),kron(speye(mesh.n[2]),av(mesh.n[1])))
  A3 = kron(speye(mesh.n[3]),kron(av(mesh.n[2]),av(mesh.n[1])))
  d1 = wMat*(A1*data[1:mesh.ne[1]])
  d2 = wMat*(A2*data[mesh.ne[1]+1:mesh.ne[1]+mesh.ne[2]])
  d3 = wMat*(A3*data[mesh.ne[1]+mesh.ne[2]+1:sum(mesh.ne)])
  dout = [d1';d2';d3']
  return dout
end

#-----------------------------------------------------------------

function CreateVTKtimeSeriesCollection(baseFname::String,mesh::TensorMesh3D,
                        t::Array{Float64,1},data...)
  #Convert mesh to vtk rectilinear format
  x,y,z = getXYZvtk(mesh)
  
  #Initialise pvd container file
  pvd = paraview_collection(baseFname)
  
  #Add each time step
  nt = length(t)
  vtk = []
  for it = 0:nt-1
    push!(vtk,vtk_grid(string(baseFname,@sprintf("_%02i",it)),x,y,z))
    for idat in data
      if (length(idat) == 2)
        if (size(idat[1],2) > 1)
          addVTKdataTensor(vtk[end],idat[1][:,it+1],idat[2],mesh)
        else
          addVTKdataTensor(vtk[end],idat[1],idat[2],mesh)
        end
      else
        if (size(idat[1],2) > 1)
          addVTKdataTensor(vtk[end],idat[1][:,it+1],idat[2],mesh,idat[3])
        else
          addVTKdataTensor(vtk[end],idat[1],idat[2],mesh,idat[3])
        end
      end
    end
    collection_add_timestep(pvd,vtk[end],t[it+1])
  end
  vtk_save(pvd)
end

function getXYZvtk(mesh::TensorMesh3D)
  x0 = mesh.x0
  h1 = mesh.h1
  h2 = mesh.h2
  h3 = mesh.h3

  x    = zeros(mesh.n[1]+1)
  x[1] = x0[1]
  y    = zeros(mesh.n[2]+1)
  y[1] = x0[2]
  z    = zeros(mesh.n[3]+1)
  z[1] = x0[3]
  for ix = 1:mesh.n[1]
    x[ix+1] = x[ix]+h1[ix]
  end
  for iy = 1:mesh.n[2]
    y[iy+1] = y[iy]+h2[iy]
  end
  for iz = 1:mesh.n[3]
    z[iz+1] = z[iz]+h3[iz]
  end
  return x,y,z
end

