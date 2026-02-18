#!/usr/bin/env python3
"""Launch the sweep viewer. Starts a local HTTP server and opens the viewer in your browser."""
import webbrowser
import http.server
import socketserver
import os

PORT = 8000
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

os.chdir(SCRIPT_DIR)

class Handler(http.server.SimpleHTTPRequestHandler):
    def end_headers(self):
        self.send_header("Cache-Control", "no-store, no-cache, must-revalidate")
        super().end_headers()

with socketserver.TCPServer(("", PORT), Handler) as httpd:
    url = f"http://localhost:{PORT}/sweep_viewer.html"
    print(f"Sweep viewer: {url}")
    webbrowser.open(url)
    print("Press Ctrl+C to stop")
    httpd.serve_forever()
