// Custom JavaScript for documentation
document.addEventListener('DOMContentLoaded', function() {
    // Initialize Mermaid diagrams
    if (typeof mermaid !== 'undefined') {
        mermaid.initialize({
            startOnLoad: true,
            theme: 'default',
            themeVariables: {
                primaryColor: '#007bff',
                primaryTextColor: '#333',
                primaryBorderColor: '#007bff',
                lineColor: '#333',
                secondaryColor: '#f8f9fa',
                tertiaryColor: '#e9ecef'
            }
        });
    }
    
    // Add copy buttons to code blocks
    document.querySelectorAll('pre code').forEach(function(block) {
        const button = document.createElement('button');
        button.className = 'copy-button';
        button.textContent = 'Copy';
        button.onclick = function() {
            navigator.clipboard.writeText(block.textContent);
            button.textContent = 'Copied!';
            setTimeout(() => button.textContent = 'Copy', 2000);
        };
        block.parentNode.appendChild(button);
    });
});

// MathJax configuration
window.MathJax = {
  tex: {
    inlineMath: [["\\(", "\\)"]],
    displayMath: [["\\[", "\\]"]],
    processEscapes: true,
    processEnvironments: true
  },
  options: {
    ignoreHtmlClass: ".*|",
    processHtmlClass: "arithmatex"
  }
};

// Initialize Mermaid
document.addEventListener('DOMContentLoaded', function() {
  if (typeof mermaid !== 'undefined') {
    mermaid.initialize({
      startOnLoad: true,
      theme: 'default'
    });
  }
});