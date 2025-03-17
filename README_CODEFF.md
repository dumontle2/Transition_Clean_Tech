# Prerequisites

- [Docker Desktop for Windows](https://docs.docker.com/desktop/setup/install/windows-install/)

# Run container in detached mode.

Note during first run it will take more time for building image

```
docker compose up -d
```

# Enter to the container

```
docker exec -it python_app bash
```

# Run the scripts

```
cd code
python estimation_weave.py
python load_data.py
python generate_policy.py
```

# Output

You can find output in the `figures` folder in your `code` folder
