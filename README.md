# AI-Opti-Solver
This is an old project from the artificial intelligence lecture. It solves a cost function with given constraints.

## Warning: unstable and not well documented
This is code is old and can be compiled but ran in a seg fault error when executing (back when i wrote this code it worked fine in those days).

## What does it do?
It reads in a taskfile where a costfunction is given (max: ... or min: ...) and try to optimize a solution vector. To do so it need constraints which have also be present in the taskfile.

### There are some defines
The `creadted_task.txt` is the result of the in-buld feature of generating a random taskfile. By setting the `#define CREATE_TASKFILE 0` to 1 you can override the random taskfile.

## To-Do
- fix issues on linux (make it runable again)
- explain here in the readme.md the calculation background
