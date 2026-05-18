# SSH into this computer, and run this from SETD_paper/
# Runs in the backround, not interupted if ssh breaks
# Output written to remote.log
nohup /usr/ds/bin/julia --project=SFS -t 12 "modelA_conv.jl" >scripts/remote.log 2>&1 </dev/null &
