# MELTING5.2.0_4Way

## Prerequisites
1. Ensure that Docker Desktop is installed and running on your machine. If not, you can download it from [here](https://www.docker.com/products/docker-desktop).

## How to Run
Follow the steps below:

1. Clone the repository:
    ```bash
   Download the DockerFile
    ```

2. Navigate into the cloned directory:
    ```bash
   docker build -t melting .
    ```

3. Run the bash script by providing the path to your specific folder which contains the input CSV file. For example:
    ```bash
   bash test.sh ./target
    ```
    **Note:** Replace `/path/to/your/folder` with the actual path to your folder containing the input CSV file.

If you encounter any issues, feel free to open a new issue in this repository.
