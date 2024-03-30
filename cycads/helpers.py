import subprocess

def which(command):
    try:
        result = subprocess.run(['which', command], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        return result.stdout.strip()
    except subprocess.CalledProcessError:
        return None


