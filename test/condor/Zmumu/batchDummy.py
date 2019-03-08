universe = vanilla
Executable = exe
should_transfer_files = YES
WhenTOTransferOutput  = ON_EXIT
notification = never
myPath = path
Output = $(myPath)/batch_$(cluster)_$(process).stdout
Error  = $(myPath)/batch_$(cluster)_$(process).stderr
Log    = $(myPath)/batch_$(cluster)_$(process).condor
Requirements = OpSys == "LINUX" && (Arch != "DUMMY" )
request_disk = 10000000
request_memory = 2100
Arguments = $(process)
Queue queue
