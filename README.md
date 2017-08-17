# THAP: A matlab Tool for HAwkes Processes

_THAP_ is a pure matlab toolbox for modeling and analysis of Hawkes process and its variants. The license is shown in [LICENSE](LICENSE).

The project was started in 2017 by Hongteng Xu and Hongyuan Zha at the School of [Computational Science and Engineering](https://www.cse.gatech.edu/) of [Georgia Institute of Technology](http://www.gatech.edu/), Atlanta, USA. 

## Quick description

The focus of _THAP_ is on learning a special kind of point processes called Hawkes process (HP). 
This toolbox is composed of multiple simulators, models, learning algorithms, and analytic and visualization tools of Hawkes processes and their variants. 
Typical functions achieved by the toolbox include:

- Simulate (multi-dimensional) Hawkes processes with smooth kernels like exponential and Gaussian kernels. 
- Segment or stitch asynchronous event sequences.
- Learn parametric or nonparametric Hawkes models by maximun likelihood estimation (MLE) and least-square estimation (LS).
- Learn typical variants of traditional Hawkes processes, e.g., mixtures of Hawkes processes (MixHP) and time-varying Hawkes processes (TVHP).
- Construct the Granger causality graph of events.
- Caculate distance between marked point processes.
- Visualize asynchronous event sequences, their intensity functions and the impact functions of Hawkes processes.
- Predict future events based on historical observations and learned models.

A set of examples can be found at 
- [https://github.com/HongtengXu/Hawkes-Process-Toolkit/](https://github.com/HongtengXu/Hawkes-Process-Toolkit/) 

The paper associated to this toolbox will come soon! 

If you use _THAP_ in publications, we would appreciate citations of related papers.

## Typical applications

The typical applications of _THAP_ include, but not limited to:

_Healthcare_: 
* Construct disease networks from patients' admission history.
* Model and predict the transitions of patients between different care units.
   
_Social science_:
* IPTV user behavior analysis: find clusters of IPTV users according to their TV program viewing records. 
* Social network analysis: model user infectivity on social networks according to users' social behaviors. 
* Talent flow modeling: model the job hopping behaviors between companies.

## Requirements and installation

_THAP_ is compatible on Windows/Linux/OSX systems and only requires MatlabR2016a or newer. 
To use _THAP_, you just need to open your Matlab and run the "setup.m" file to add paths of necessary functions.


## References and citation

According to the functions used in your work, you may want to cite the following papers:

_Simulation methods of Hawkes processes_:
* Ogata, Yosihiko. "[On Lewis' simulation method for point processes.](http://ieeexplore.ieee.org/abstract/document/1056305/)" IEEE Transactions on Information Theory 27.1 (1981): 23-31. 
* Møller, Jesper, and Jakob G. Rasmussen. "[Approximate simulation of Hawkes processes.](https://link.springer.com/article/10.1007%2Fs11009-006-7288-z?LI=true)" Methodology and Computing in Applied Probability 8.1 (2006): 53-64. 
* Dassios, Angelos, and Hongbiao Zhao. "[Exact simulation of Hawkes process with exponentially decaying intensity.](http://emis.ams.org/journals/EJP-ECP/article/download/2717/2717-14258-1-PB.pdf)" Electronic Communications in Probability 18 (2013). 

_Learning algorithms of Hawkes processes_:
* Zhou, Ke, Hongyuan Zha, and Le Song. "[Learning social infectivity in sparse low-rank networks using multi-dimensional hawkes processes.](http://proceedings.mlr.press/v31/zhou13a.html)" Artificial Intelligence and Statistics. 2013. 
* Zhou, Ke, Hongyuan Zha, and Le Song. "[Learning triggering kernels for multi-dimensional hawkes processes.](http://proceedings.mlr.press/v28/zhou13.pdf)" Proceedings of the 30th International Conference on Machine Learning (ICML-13). 2013. 
* Luo, Dixin, et al. "[Multi-Task Multi-Dimensional Hawkes Processes for Modeling Event Sequences.](http://www.aaai.org/ocs/index.php/IJCAI/IJCAI15/paper/download/11299/11253)" IJCAI. 2015. 
* Eichler, Michael, Rainer Dahlhaus, and Johannes Dueck. "[Graphical modeling for multivariate hawkes processes with nonparametric link functions.](http://onlinelibrary.wiley.com/doi/10.1111/jtsa.12213/full)" Journal of Time Series Analysis 38.2 (2017): 225-242. 

_Time-varying Hawkes processes_: 
* Xu, Hongteng, Dixin Luo, and Hongyuan Zha. "[Learning Hawkes Processes from Short Doubly-Censored Event Sequences.](https://arxiv.org/pdf/1702.07013.pdf)" arXiv preprint arXiv:1702.07013 (2017). 
* Xu, Hongteng, et al. "[Patient flow prediction via discriminative learning of mutually-correcting processes.](http://ieeexplore.ieee.org/abstract/document/7593326/)" IEEE transactions on Knowledge and Data Engineering 29.1 (2017): 157-171. 

_Granger causality analysis of event sequences_: 
* Xu, Hongteng, Mehrdad Farajtabar, and Hongyuan Zha. "[Learning granger causality for hawkes processes.](http://proceedings.mlr.press/v48/xuc16.pdf)" International Conference on Machine Learning. 2016. 

_Clustering analysis of event sequences_:
* Xu, Hongteng, and Hongyuan Zha. "[A Dirichlet Mixture Model of Hawkes Processes for Event Sequence Clustering.](https://arxiv.org/pdf/1701.09177.pdf)" arXiv preprint arXiv:1701.09177 (2017).
* Iwayama, Koji, Yoshito Hirata, and Kazuyuki Aihara. "[Definition of distance for nonlinear time series analysis of marked point process data.](http://www.sciencedirect.com/science/article/pii/S037596011631564X)" Physics Letters A 381.4 (2017): 257-262.

_The usage of data_:
* IPTV data: Luo, Dixin, et al. "[You are what you watch and when you watch: Inferring household structures from iptv viewing data.](http://ieeexplore.ieee.org/abstract/document/6717182/)" IEEE Transactions on Broadcasting 60.1 (2014): 61-72. 
* Linkedin data: Xu, Hongteng, Dixin Luo, and Hongyuan Zha. "[Learning Hawkes Processes from Short Doubly-Censored Event Sequences.](https://arxiv.org/pdf/1702.07013.pdf)" arXiv preprint arXiv:1702.07013 (2017). 
