FROM python:2.7.17

# Set environment variables to avoid prompts
ENV DEBIAN_FRONTEND=noninteractive

# Update package lists and install dependencies
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    python-pip \
    python-dev \
    libatlas-base-dev \
    build-essential \
    && rm -rf /var/lib/apt/lists/*

# Upgrade pip
RUN pip install --upgrade pip

# Install Python packages
RUN pip install numpy
RUN pip install scipy
RUN pip install scipy-weave
RUN pip install matplotlib
RUN pip install seaborn

# Set the working directory
WORKDIR /app

# Default command
CMD ["python"]
