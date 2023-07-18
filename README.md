# MELTING5.2.0_4Way

## Prerequisites
1. Ensure that Docker Desktop is installed and running on your machine. If not, you can download it from [here](https://www.docker.com/products/docker-desktop).

## How to Run
Follow the steps below:

1. 
    ```bash
   Download the Dockerfile and test.sh files
    ```

2. Run the below in terminal after going to the directory where the downloaded dockerfile is(The name of Dockerfile should be Dockerfile and also remove the .txt extension):
    ```bash
   docker build -t melting .
    ```

4. Run the bash script by providing the path to your specific folder which contains the input CSV file. For example:
    ```bash
   bash test.sh ./target
    ```
    **Note:** Replace `/path/to/your/folder` with the actual path to your folder containing the input CSV file.
   
5. You will find the output.csv in the ./target folder
If you encounter any issues, feel free to open a new issue in this repository.
