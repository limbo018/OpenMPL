# SimpleMPL
---------

**SimpleMPL** stands for simple multiple patterning lithography framework

---------
## Features
 * Contact or metal layer decomposition 
 * No stitching
 * Support 3 or 4 coloring 
 * Density control
 * Multi-threading
 * Small memory usage
 * Multiple algorithms: 
     ILP (Gurobi or Lemon CBC), 
     SDP (Csdp API), 
     LP  (Gurobi API)

The Csdp API used in SimpleMPL has been modified and built for threading safety at high level. 

---------
## Authors

|  Name              | Affiliation                |  email                            |
| ------------------ | -------------------------- | --------------------------------- |
| Yibo Lin           | ECE Dept, UT Austin        | yibolin@cerc.utexas.edu           |
| Bei Yu             | CSE Dept, CUHK             | byu@cse.cuhk.edu.hk               |
| David Z. Pan       | ECE Dept, UT Austin        | dpan@ece.utexas.edu               |
