# THAP: A matlab Tool for HAwkes Processes

_THAP_ is a pure matlab toolbox for modeling and analysis of Hawkes process and its variants. The license is shown in [LICENSE](LICENSE).

The project was started in 2017 by Hongteng Xu and Hongyuan Zha at the School of Computational Science and Engineering (https://www.cse.gatech.edu/) of [Georgia Institute of Technology](http://www.gatech.edu/), Atlanta, USA. 

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

## Requirements and Installation

_THAP_ is compatible on Windows/Linux/OSX systems and only requires MatlabR2016a or newer. 
To use _THAP_, you just need to open your Matlab and run the "setup.m" file to add paths of necessary functions.


## References and Citation

According to the functions you used in your work, you may want to cite the following papers:

_Learning algorithms of Hawkes processes_:
* Zhou, Ke, Hongyuan Zha, and Le Song. "Learning social infectivity in sparse low-rank networks using multi-dimensional hawkes processes." Artificial Intelligence and Statistics. 2013. [http://proceedings.mlr.press/v31/zhou13a.html](http://proceedings.mlr.press/v31/zhou13a.html)
* Zhou, Ke, Hongyuan Zha, and Le Song. "Learning triggering kernels for multi-dimensional hawkes processes." Proceedings of the 30th International Conference on Machine Learning (ICML-13). 2013. [http://proceedings.mlr.press/v28/zhou13.pdf](http://proceedings.mlr.press/v28/zhou13.pdf)
* Luo, Dixin, et al. "Multi-Task Multi-Dimensional Hawkes Processes for Modeling Event Sequences." IJCAI. 2015.[http://www.aaai.org/ocs/index.php/IJCAI/IJCAI15/paper/download/11299/11253](http://www.aaai.org/ocs/index.php/IJCAI/IJCAI15/paper/download/11299/11253)
* Eichler, Michael, Rainer Dahlhaus, and Johannes Dueck. "Graphical modeling for multivariate hawkes processes with nonparametric link functions." Journal of Time Series Analysis 38.2 (2017): 225-242. [http://onlinelibrary.wiley.com/doi/10.1111/jtsa.12213/full](http://onlinelibrary.wiley.com/doi/10.1111/jtsa.12213/full)

_Time-varying Hawkes processes_: 
* Xu, Hongteng, Dixin Luo, and Hongyuan Zha. "Learning Hawkes Processes from Short Doubly-Censored Event Sequences." arXiv preprint arXiv:1702.07013 (2017). [https://arxiv.org/pdf/1702.07013.pdf](https://arxiv.org/pdf/1702.07013.pdf)
* Xu, Hongteng, et al. "Patient flow prediction via discriminative learning of mutually-correcting processes." IEEE transactions on Knowledge and Data Engineering 29.1 (2017): 157-171. [http://ieeexplore.ieee.org/abstract/document/7593326/](http://ieeexplore.ieee.org/abstract/document/7593326/)

_Granger causality analysis of event sequences_: 
* Xu, Hongteng, Mehrdad Farajtabar, and Hongyuan Zha. "Learning granger causality for hawkes processes." International Conference on Machine Learning. 2016. [http://proceedings.mlr.press/v48/xuc16.pdf](http://proceedings.mlr.press/v48/xuc16.pdf)

_Clustering analysis of event sequences_:
* Xu, Hongteng, and Hongyuan Zha. "A Dirichlet Mixture Model of Hawkes Processes for Event Sequence Clustering." arXiv preprint arXiv:1701.09177 (2017). [https://arxiv.org/pdf/1701.09177.pdf](https://arxiv.org/pdf/1701.09177.pdf)
* Iwayama, Koji, Yoshito Hirata, and Kazuyuki Aihara. "Definition of distance for nonlinear time series analysis of marked point process data." Physics Letters A 381.4 (2017): 257-262. [http://www.sciencedirect.com/science/article/pii/S037596011631564X](http://www.sciencedirect.com/science/article/pii/S037596011631564X)
