 #!/bin/bash
   #SBATCH   --partition=normal            # submit   to the normal(default) partition
   #SBATCH   --job-name=t-test             # name the job
   #SBATCH   --output=t-test-%j.out        # write stdout/stderr   to named file
   #SBATCH   --error=t-test-%j.err      
   #SBATCH   --time=3-00:00:00
   #SBATCH   --nodes=1                     # Request N nodes
   #SBATCH   --cpus-per-task=48            # Request n   cores per node
   #SBATCH   --mem=192GB             # Request nGB RAM per core

   ./main.out