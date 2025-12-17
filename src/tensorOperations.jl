export 
  axis_core_tensor, 
  add_slice_all!, 
  multiply_slice_all!, 
  swap_slice_all!, 
  matrix_tensor_multiply!, 
  matrix_tensor_multiply, 
  matrix_tensor_congruence!, 
  matrix_tensor_congruence

function axis_core_tensor(_d)
  C = zeros(QQ,_d,_d,_d);
  for al in (1:_d)
    for be in (1:_d)
      for ga in (1:_d)
        if al == be && be == ga
          C[al,be,ga] = QQ(1,1)
        end
        if (al < be && be == ga)||(al == be && be < ga)
          C[al,be,ga] = QQ(3,1)
        end
        if al < be && be < ga
          C[al,be,ga] = QQ(6,1)
        end
      end
    end
  end # QQ.(C) == tensor_sequence(sig_axis(TTSd))[4];
  return C
end

function add_slice_all!(G,a,i::Int,j::Int)
  G[j,:,:] += a.*G[i,:,:]
  G[:,j,:] += a.*G[:,i,:]
  G[:,:,j] += a.*G[:,:,i]
  return G
end

function multiply_slice_all!(G,a,i::Int)
  G[i,:,:] = a.*G[i,:,:]
  G[:,i,:] = a.*G[:,i,:]
  G[:,:,i] = a.*G[:,:,i]
  return G
end

function swap_slice_all!(G,i::Int,j::Int)
  temp = G[i,:,:]
  G[i,:,:] = G[j,:,:]
  G[j,:,:] = temp
  temp = G[:,i,:]
  G[:,i,:] = G[:,j,:]
  G[:,j,:] = temp
  temp = G[:,:,i]
  G[:,:,i] = G[:,:,j]
  G[:,:,j] = temp
  return G
end

function matrix_tensor_multiply!(matrix::Array, tensor::Array{T,3}, index::Int) where T
    @assert 1 ≤ index ≤ 3 "Index must be 1, 2, or 3."
    @assert size(tensor, index) == size(matrix, 2) "Dimension mismatch."
    if index == 1
        # (n, a, b) → (m, a, b)
        resh = reshape(tensor, size(tensor,1), :)
        result = matrix * resh
        copyto!(tensor, reshape(result, size(matrix,1), size(tensor,2), size(tensor,3)))
    elseif index == 2
        # (a, n, b) → (a, m, b)
        # permute to (n, a, b)
        tensor_perm = permutedims(tensor, (2,1,3))
        resh = reshape(tensor_perm, size(tensor,2), :)
        result = matrix * resh
        result = reshape(result, size(matrix,1), size(tensor,1), size(tensor,3))
        copyto!(tensor, permutedims(result, (2,1,3)))
    else  # index == 3
        # (a, b, n) → (a, b, m)
        tensor_perm = permutedims(tensor, (3,1,2))
        resh = reshape(tensor_perm, size(tensor,3), :)
        result = matrix * resh
        result = reshape(result, size(matrix,1), size(tensor,1), size(tensor,2))
        copyto!(tensor, permutedims(result, (2,3,1)))
    end
    return tensor
end

function matrix_tensor_congruence!(matrix::Array, tensor::Array{T,3}) where T
    @assert size(tensor,1) == size(matrix,2)
    @assert size(tensor,2) == size(matrix,2)
    @assert size(tensor,3) == size(matrix,2)
    matrix_tensor_multiply!(matrix, tensor, 1)
    matrix_tensor_multiply!(matrix, tensor, 2)
    matrix_tensor_multiply!(matrix, tensor, 3)
    return tensor
end

function matrix_tensor_multiply(matrix::Array, tensor::Array, index::Int)
    """
    Perform matrix-tensor multiplication along the specified axis of a k-tensor.
    Parameters:
        matrix::Array: The matrix to multiply (2D array).
        tensor::Array: The input tensor (k-dimensional array).
        index::Int: The axis of the tensor to perform the multiplication.
    Returns:
        A new tensor resulting from the multiplication.
    """
    # Ensure the dimensions match for contraction
    @assert 1 ≤ index ≤ ndims(tensor) "Index must be a valid axis of the tensor."
    @assert size(tensor, index) == size(matrix, 2) "Matrix and tensor dimensions do not align!"
    # Move the target axis to the first dimension
    perm = [index; 1:index-1; index+1:ndims(tensor)]
    tensor_perm = permutedims(tensor, perm)
    # Reshape the tensor to (size along index, rest)
    reshaped_tensor = reshape(tensor_perm, size(tensor, index), :)
    # Perform matrix multiplication
    result = matrix * reshaped_tensor  # Result shape: (matrix rows, remaining dimensions)
    # Reshape back to the original dimensions
    result_shape = (size(matrix, 1), size(tensor)[[1:index-1; index+1:end]]...)
    result = reshape(result, result_shape)
    # Undo the permutation to restore original axes order
    inv_perm = invperm(perm)
    return permutedims(result, inv_perm)
end

function matrix_tensor_congruence(matrix::Array, tensor::Array)
    k = length(size(tensor))
    if k == 0 
      return tensor
    end 
    res = tensor
    #res = Array{eltype(tensor)}(undef, size(matrix, 1)*ones(Int,k)...)
    for i in (1:k)
        @assert size(tensor, i) == size(matrix, 2) "Matrix and tensor dimensions do not align!"
        res = matrix_tensor_multiply(matrix, res, i)
    end
    return res
end

