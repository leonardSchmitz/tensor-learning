# optimization-of-tensor-learning

This is the companion package to the article  "Optimization of Tensor Learning" [[1]](#1) 
in the computer algebra system [OSCAR](https://www.oscar-system.org). 

## Reference

<a id="1">[1]</a>
C. Amendola, L. Schmitz </br>
[Learning barycenters from signature matrices](https://arxiv.org/abs/2509.07815v1) </br>
arXiv:2509.07815, 2025. 

## Contact

**Authors**

- **Leonard Schmitz**  
  Technische Universität Berlin, Germany 
  Email: [lschmitz[at]math.tu-berlin.de](mailto:lschmitz@math.tu-berlin.de)

## Set up
1. Clone the repository by running
   ```
   git clone https://github.com/leonardSchmitz/optimization-of-tensor-learning.git
   ```
2. Install [OSCAR](https://www.oscar-system.org)
3. Open julia in the project folder via
   ```
   julia 
   ``` 
4. The first time using the package, use
   ```
   using Pkg 
   Pkg.instantiate()
   ```
3. Start the package with  
   ```
   using Oscar 
   using optTensLearning
   ```

## File dictionary 
We provide the source code for several compuations presented in [[1]](#1).
The folloing list helps the reader to find the relevant source file.  

| Reference in [[1]](#1)      | name of file                         |
|-----------------------------|--------------------------------------|
| Ex 2.2,2.4,3.3,3.5,4.3,4.5  | `examples/running_example.jl`        |
| Proposition 4.1             | `proofs/Proposition41.jl`            |
| Proposition 4.4             | `proofs/Proposition44.jl`            |
| Table 1 (ours)              | `tests/learning_tests.jl`            |
| Table 1 (GBF4)              | `tests/GB_learning_tests.jl`         |


## Efficiency remarks

The purpose of this package is to support research, verify compuitational results, and enable experimentation. 
Several solutions have been kept simple and are not optimal from an implementation perspective. 
Achieving optimality is the next step. 
We report performance timings on a MacBook Air (1.6 GHz Dual-Core Intel Core i5), see Table 1 in [[1]](#1). 



