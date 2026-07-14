build with

```bash/fish/zsh/pwsh
gcc main.c -O0 -g -Og -march=native -fsanitize=address,leak,undefined -Wall -Wextra -Wpedantic -Werror -o main
```
