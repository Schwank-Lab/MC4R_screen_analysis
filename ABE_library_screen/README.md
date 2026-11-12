## Setup Conda Environment

1. **Create environment (first time):**
	```bash
	conda create -n mc4r_screen python=3.10
	conda activate mc4r_screen
	```

2. **Export environment (for backup or sharing):**
	```bash
	conda env export > MC4R_screen_env.yml
	```

3. **Recreate environment from file:**
	```bash
	conda env create -n mc4r_screen -f MC4R_screen_env.yml
	conda activate mc4r_screen
	```

## How to Run

After activating the environment, run the pipeline:
```bash
python run_pipeline.py
```