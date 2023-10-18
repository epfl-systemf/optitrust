
OptiTrust is a tool for user-guided, source-to-source transformations.

It is a research prototype, in active development. We hope to reach by the end
of 2023 a stage of robustness and documentation that will make is feasible
for external users to play with OptiTrust and contribute new transformations.

In the meanwhile, if you are interested in a demo, please get in touch with @charguer.


# OptiTrust Installation


--------------------------------------------------------------------------------
## Installation

# Install OCaml and system packages

It takes about 30 minutes to install the required OCaml software.

Installation of system packages:

```sh
   sudo apt-get install opam clang-format libclang-dev llvm-dev libomp-dev pkg-config zlib1g-dev
   # optional, only if you prefer using `meld` over `code -d`:
   sudo apt-get install meld
```

Installation of OCaml ecosystem:

```sh
   opam init
   opam switch create 4.14.1
   opam pin add menhirLib 20210419
   opam pin add pprint 20220103
   opam pin add clangml 4.8.0
   opam install dune refl clangml pprint menhir menhirLib base64 ocamlbuild ocaml-lsp-server
   # next line used only for generating the documentation of OptiTrust:
   opam install odoc lambdasoup
   # then in any case execute the last line below
   eval $(opam env)
```

# Install libraries for parsing

Then, you need to either export the environment variable OPTITRUST by
executing "export OPTITRUST=`pwd`", from the optitrust folder,
or more conveniently execute the command:

```sh
  sudo make install_compcert_stdlib"
```

which essentially performs a `sudo install src/c/compcert_parser/include/*.h usr/local/lib/compcert`

# Test your installation from the command line

Checking your installation of OptiTrust is working:
```sh
   make tests
```

# Configure VScode (or VSCodium) for interactive usage

You can install either VSCode or VSCodium (more open).

Installation of VSCode (see https://code.visualstudio.com/download for alternatives)
```sh
  sudo snap install --classic code
```

Installation of VSCodium:

```sh
  sudo snap install codium --classic.
```

Alternative without use of snap:

```sh
   wget -qO - https://gitlab.com/paulcarroty/vscodium-deb-rpm-repo/raw/master/pub.gpg \
    | gpg --dearmor | sudo dd of=/usr/share/keyrings/vscodium-archive-keyring.gpg
   echo 'deb [ signed-by=/usr/share/keyrings/vscodium-archive-keyring.gpg ] https://download.vscodium.com/debs vscodium main' \
    | sudo tee /etc/apt/sources.list.d/vscodium.list
   sudo apt update
   sudo apt install codium
```

In the rest of this tutorial, we assume VSCode. If you are using VSCodium,
you probably have an existing alias from 'code' to 'codium'. If not, you
may want to add into your `~/.bashrc` the line:
```bash
   alias code='codium'
```
(or use `sudo ln -s /usr/bin/codium /usr/bin/code`).

Recommended: installation of Chromium browser, which is very fast for
opening up the "diff" pages generated by OptiTrust:

```sh
   sudo snap install chromium
```

--------------------------------------------------------------------------------
## Configuration of VSCode

The point is to set up keybindings for interactive development of tranformation
scripts. For example, F6 shows the diff associated with the transformation
covering the cursor position in the VSCode editor.

### Open VSCode

From the folder that contains present README file, open VSCode using:

```sh

   code . &
   # or
   codium . &

```

### Install the OCaml platform:

Open the extension pannel (`ctrl+shift+x`) and look for the "Ocaml platform" extension.
Alternatively, use the quick open prompt (`ctrl+p`), then paste `ext install ocamllabs.ocaml-platform`.

### Install the OptiTrust shortcuts for VSCode

In VSCode, open the file `~/.config/Code/User/keybindings.json`.
For VSCodium, this file is located at `~/.config/VSCodium/User/keybindings.json`.
If you have an empty file, paste the following contents.
If you have a nonempty file, copy the inner contents of the outermost braces,
and merge that contents just before the final closing brace of the existing file.


```jsonc
  // OptiTrust keybindings
  {
    "key": "f5",
    "command": "workbench.action.tasks.runTask",
    "args": "Redo last view command",
    "when": "config.optitrust.enableKeybindings"
  },
  {
    "key": "shift+f5",
    "command": "workbench.action.tasks.runTask",
    "args": "View trace",
    "when": "config.optitrust.enableKeybindings"
  },
  {
    "key": "f6",
    "command": "workbench.action.tasks.runTask",
    "args": "View diff",
    "when": "config.optitrust.enableKeybindings"
  },
  {
    "key": "alt+f6",
    "command": "workbench.action.tasks.runTask",
    "args": "View big step diff",
    "when": "config.optitrust.enableKeybindings"
  },
  {
    "key":"ctrl+shift+f6",
    "command": "workbench.action.tasks.runTask",
    "args": "View diff for ast encoding",
    "when": "config.optitrust.enableKeybindings"
  },
  // For working with unit tests
  {
    "key": "f10",
    "command": "workbench.action.tasks.runTask",
    "args": "Rerun the last-tried test(s)",
    "when": "config.optitrust.enableKeybindings"
  },
  {
    "key": "ctrl+f10",
    "command": "workbench.action.tasks.runTask",
    "args": "Run the current test",
    "when": "config.optitrust.enableKeybindings"
  },
  {
    "key": "shift+f10",
    "command": "workbench.action.tasks.runTask",
    "args": "Run all the tests",
    "when": "config.optitrust.enableKeybindings"
  },
  {
    "key": "alt+f10",
    "command": "workbench.action.tasks.runTask",
    "args": "Open unit test ML and CPP files",
    "when": "config.optitrust.enableKeybindings"
  },
  // For working with long transformation scripts
  {
    "key": "f7",
    "command": "workbench.action.tasks.runTask",
    "args": "View diff from intermediate state",
    "when": "config.optitrust.enableKeybindings"
  },
  {
    "key": "alt+f7",
    "command": "workbench.action.tasks.runTask",
    "args": "View big step diff from intermediate state",
    "when": "config.optitrust.enableKeybindings"
  },
  {
    "key": "ctrl+f7",
    "command": "workbench.action.tasks.runTask",
    "args": "Save intermediate state",
    "when": "config.optitrust.enableKeybindings"
  },
  // For killing a task, type 'ctrl+k' twice, then 'enter'
  {
     "key": "ctrl+k ctrl+k",
     "command": "workbench.action.tasks.terminate",
     "when": "config.optitrust.enableKeybindings"
  },

```

Note: the shortcuts refer to tasks that are defined in `.vscode/tasks.json`,
which can be opened by pressing `ctrl+p` then typing `tasks.json`.

You can then adjust the keybindings to match your personal preferences,
by changing the key parameter.


Note: the intermediate state shortcuts (`F7`) are useful for saving checkpoints.
Press `ctrl+F7` to compute the result of a transformation at a given line
of the script. Then, type `F7` on a line further in the file.
The result should be the same as with typing `F6`, except that it runs
faster because the execution proceeds not from the beginning of the script
but instead from the line of the snapshot taken using `ctrl+F7`.

### Deactivate conficting Ubuntu binding

IMPORTANT: on Ubuntu, `Alt+F6`, `Alt+F7` etc. are bound to window manipulation operations,
e.g. resize. You can either modify the shortcut, or (easier) deactivate the Ubuntu binding.
To that end you may use the following commands:

```sh
  sudo apt-get install dconf-editor
  gsettings set org.gnome.desktop.wm.keybindings begin-move []
  gsettings set org.gnome.desktop.wm.keybindings begin-resize []
  gsettings set org.gnome.desktop.wm.keybindings cycle-group []
  gsettings set org.gnome.desktop.wm.keybindings cycle-group-backward []
  gsettings set org.gnome.desktop.wm.keybindings toggle-maximized []
```

Alternatively, you can open the settings panel, the keyboard menu, the
shortcut submenu, then in the search bar type `Alt+F` (or go to the
'window' group), then select an action, and type `Backspace` to disable
the shortcut, then click on the `save` button.


--------------------------------------------------------------------------------
## Using OptiTrust in VSCode


### Test your installation

We are at last ready to test the installation on a case study.

In VScode, open `case_studies/matmul/matmul.ml` (using e.g. `ctrl+p`).
Place your cursor on the line starting with `!! List.iter tile`.
Type `F6`.
You should see a diff opening up in a browser.

Another shortcut to try is `shift+F5`, on any line of the `matmul.ml` file
After a dozen seconds, you should see the full transformation trace for the
matrix-multiply case study.

### Troubleshooting

If you don't see a diff, possible issues include:
   - Your shortcut is not set up correctly; when the shortcut is pressed,
   a terminal should open to show the output of the task.
     If needed, to investigate whether key bindings are set up properly,
     in VScode type `ctrl+K` immediately followed by `ctrl+s`, then type `alt+k`,
     then type `F6` and see whether you see "Tasks: run tasks" a entry.
   - the compilation failed due to incorrect setup; you should see error
     messages in the terminal.



--------------------------------------------------------------------------------
## Optional suggestions

In your `~/.bashrc` you can add an alias to make it easier to invoke the tester.

```sh
  alias t='~/optitrust/tester'
```


--------------------------------------------------------------------------------
## Optional tools

See `INSTALL_EXTRA.md` for a list of additional useful tools for program optimization.

See `VSCODE_CUSTOMIZE.md` for useful tips for using VScode or VScodium.



