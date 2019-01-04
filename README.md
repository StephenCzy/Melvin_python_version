SimpleMelvinExamples from Dr. Mario Krenn 

https://mariokrenn.wordpress.com/

SimpleHOMExample.nb

This program shows how to work with quantum states, and how the symbolic transformations work.

CalcSRV.nb

This is a full version which searches for 3-particle high-dimensionally entanged states with existing optical elements.

More Information
Mario Krenn, Mehul Malik, Robert Fickler, Radek Lapkiewicz, and Anton Zeilinger. Automated search for new quantum experiments. Physical Review Letters, 116(9):090405, 2016.

BibTex:
@article{krenn2016automated,
  title={Automated search for new quantum experiments},
  author={Krenn, Mario and Malik, Mehul and Fickler, Robert and Lapkiewicz, Radek and Zeilinger, Anton},
  journal={Physical review letters},
  volume={116},
  number={9},
  pages={090405},
  year={2016},
  publisher={APS}
}


# Quantum Information for Developers
![Docker Build Status](https://img.shields.io/docker/build/jrottenberg/ffmpeg.svg)

# Installation
```
	docker pull westernmagic/quid:2018
	# See docker create --help for details
	docker create \
		--interactive \
		--tty \
		--publish 8888:8888 \
		--restart unless-stopped \
		--volume $(pwd):/root/share \
		--name quid_2018 \
		westernmagic/quid:2018
```

# Running
```
	docker start  quid_2018
	docker attach quid_2018
```
