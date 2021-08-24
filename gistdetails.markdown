---
layout: page
title: List of PRs
permalink: /gist/
---

# GSoC '21 with GeomScale

## Project: Monte Carlo Integration
Integration is a fundamental problem in mathematics and computer science with many applications that span the 
whole spectrum of sciences and engineering. It appears, for example, in problems in statistics, biology, and 
economics, to name a few concrete application areas. Integration of functions non-rectangular domains like a 
general convex body or say a polytope turns our to be a difficult task to solve directly. 
Also, for higher dimensional integral calculation around defined integration limits can be more harder 
task while dealing with so many variables. To rule out this difficulty, Monte Carlo Integration is used. 
The overview can be found [here](https://en.wikipedia.org/wiki/Monte_Carlo_integration#Overview) for this algorithm.

Works done in Pull requests mentioned

### 1. **Simple Monte Carlo Integration**<br>
This PR is to add two functions in volesti's library to integrate functions around polytopes in N-dimensional space. 
<br>
#### **[Link to the blog post](https://surajchoubey.github.io/gsoc21/simple-mc-integration/)**<br>
#### **[Link to the PR](https://github.com/GeomScale/volume_approximation/pull/163)**<br>
#### Status: **Open + Approved**<br>
<br>
### 2. **Lovasz Vempala Integration**<BR>
This PR is to add a single function based on Lovasz Vempala research paper (Pg-7) to integrate logconcave functions around polytope in the N-dimensional space. It has substantial amount of more accuracy than Simple Monte Carlo Integration as mentioned.
<br>
#### **[Link to the blog post](https://surajchoubey.github.io/gsoc21/lv-mc-integration/)**<br>
#### **[Link to the PR](https://github.com/GeomScale/volume_approximation/pull/170)**<br>
#### Status: **Open**<br>