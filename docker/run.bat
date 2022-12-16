docker run 	--rm -it 	^
		-v %CD%\work:/home/NOMAD/work 	^
		-v %CD%/..:/NOMAD		^
		nomad_docker %*
