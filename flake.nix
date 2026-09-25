{
  description = "Two one-parameter families of rank-17 elliptic K3 surfaces (PARI/GP)";
  inputs.nixpkgs.url = "github:NixOS/nixpkgs/nixpkgs-unstable";
  outputs = { self, nixpkgs }:
    let
      forAll = f: nixpkgs.lib.genAttrs [ "x86_64-linux" "aarch64-linux" "x86_64-darwin" "aarch64-darwin" ]
        (s: f nixpkgs.legacyPackages.${s});
    in {
      devShells = forAll (pkgs: { default = pkgs.mkShell { packages = [ pkgs.pari ]; }; });
    };
}
