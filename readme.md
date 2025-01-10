# LB-PEAKS

"LB-PEAKS.py": This is the official implementation of the LB-PEAKS scheme in Python programming language. 

## Option

- [/n|-n|n]: Specify that the following option is the value of n (default: 256). 

- [/m|-m|m]: Specify that the following option is the value of m (default: 9728). 

- [/q|-q|q]: Specify that the following option is the value of q (default: 4093). 

- [/h|-h|h|/help|--help|help]: Show this help information. 

## Format

- python "LB-PEAKS.py" [/n|-n|n] n [/m|-m|m] m [/q|-q|q] q

- python "LB-PEAKS.py" [/h|-h|h|/help|--help|help]

## Example

- python "LB-PEAKS.py"

- python "LB-PEAKS.py" /n 256 /m 9728 /q 4093

- python "LB-PEAKS.py" --help

## Exit code
- 0: The Python script finished successfully. 
- 1: The Python script finished not passing all the verifications. 
- -1: The Python script received unrecognized commandline options. 

## Note

1) All the commandline options are optional and not case-sensitive. 

2) The commandline parameters will be appended to the parameter list specified by the user within the script. 

3) The parameters n, m, and q should be positive integers that are greater than 1. 

4) The parameters n and m should meet the requirement that "2n | m". Otherwise, they will be set to their default values respectively. 

## Citation

This paper is named ``Lattice-based Public Key Encryption with Aorized Keyword Search: Construction, Implementation, and Applications`` at first. 

Subsequently, it is named ``Post-Quantum Searchable Encryption Supporting User-Authorization for Outsourced Data Management``. 

Please use the following BibTeX for citations. 

```
@inproceedings{xu2024post,
  title={Post-Quantum Searchable Encryption Supporting User-Authorization for Outsourced Data Management},
  author={Xu, Shiyuan and Cao, Yibo and Chen, Xue and Guo, Yu and Yang, Yuer and Guo, Fangda and Yiu, Siu-Ming},
  booktitle={Proceedings of the 33rd ACM International Conference on Information and Knowledge Management},
  pages={2702--2711},
  year={2024}
}
```
